# Non-cyclotomic ring AKS variant — results

**Construction (S37 / Run 37, FOCUS-1 sub-attack 2 from `TODO.md`)**.
Builds and probes the AKS congruence

```
   (x + b)^n  ≡  x^n + b   (mod  f(x), n)
```

with `f(x) = x^d + a` for small `d` and small `a`, instead of the
cyclotomic `f(x) = x^r - 1` used by standard AKS. S47 (CLOSED_PATHS line
690) closed the cyclotomic-CRT shortcut: the AKS-prescribed `r` is
empirically prime, so `x^r - 1 = (x-1)·Phi_r(x)` and the non-trivial
factor still has dimension `r-1`.

The aim of this experiment was to see whether a NON-cyclotomic
modulus (Eisenstein-flavoured, `x^d + a` with `a ≥ 2`) might give a
deterministic primality test at much smaller polynomial degree, and
whether the irreducibility certificate could itself live in TC^0.

## Verdict

**PARTIAL — interesting empirical finding, but the construction does NOT
bypass the growing-dim MPOW barrier and does NOT give a proof of
correctness.**

Failure-mode tags: **(E) Equivalence + (C) Circularity**.

* (E) Working in `Z_n[x]/(x^d + a)` with `d = polylog(n)` reduces to
  matrix powering of `d × d` matrices, the *same* primitive that
  standard AKS reduces to. The Mereghetti-Palano TC^0 result for
  fixed-`k` MPOW does not extend, and the **growing-dim MPOW barrier
  (E5.3, the only open Chain-E frontier)** is unaffected.
* (C) Irreducibility certification of `x^d + a` over `Z_n` requires a
  prime factor `p | n` (Capelli's theorem applies over `F_p`, not
  `Z_n`). Producing such `p` is at least as hard as primality testing —
  circular.

## What did get found

A clean structural separation between cyclotomic and non-cyclotomic
choices on Carmichael numbers — see Q6 below. This sharpens
CLOSED_PATHS line 690 and may have small independent expository value.

## Empirical results

### Q1/Q2 — false-positive rate at small d on `n ∈ [4, 200]`

```
  Aggregate sweep (b=1, a=1):
    d=2:  primes 44/44   composites  0/153   FP = 0.0000
    d=3:  primes 44/44   composites  4/153   FP = 0.0261
    d=4:  primes 44/44   composites  0/153   FP = 0.0000
    d=5:  primes 44/44   composites  2/153   FP = 0.0131
    d=6:  primes 44/44   composites  0/153   FP = 0.0000
    d=8:  primes 44/44   composites  0/153   FP = 0.0000

  Aggregate sweep (b ∈ {1,2,3,4} all must pass, a=1):
    d=2..8: composite false-positives = 0/153 in every case.
```

So at small `n`, the test correctly catches every composite once a
handful of `b` values are tried. But `n ≤ 200` is below the Carmichael
threshold (561), so this is not yet adversarial.

### Q3 — varying `a` at fixed small `d`, `n ∈ [4, 300]`

```
  d=3,  b=1:
    a=1: FP = 0.0211   a=2: FP = 0.0000   a=3..6: FP ∈ {0.017..0.030}
  d=4,  b=1:
    a=1: FP = 0.0000   a=2: FP = 0.0000   a=5: FP = 0.0211  (others 0)
```

`a = 2` consistently gives the cleanest separation. `a = 1` is the
"cyclotomic-suspect" case (`x^d + 1` is the cyclotomic polynomial
`Phi_{2d}` for `d` a power of two; otherwise it factors over Q).

### Q4 — irreducibility cost (Capelli's theorem)

For `x^d + a` over `F_p`, irreducibility iff every prime `q | d` admits
"`-a` is not a `q`-th power in `F_p^*`", plus (if `4 | d`) "`-4a` not a
4th power". Each test is one modular exponentiation, in TC^0
(Hesse-Allender-Barrington 2002). Number of tests = `omega(d) =
O(log d / log log d) = O(log log n)` for `d = polylog(n)`.

**=> Capelli IS in TC^0 ASSUMING `p` is supplied.** But for primality
testing on input `n`, `p | n` is unknown — providing it is the
hard part of the problem. Closure mode (C) Circularity.

Empirical irreducibility frequency over random small primes:

```
  d=  2  irreducible fraction = 0.530
  d=  3                       = 0.288
  d=  4                       = 0.530
  d=  5                       = 0.179
  d=  6                       = 0.144
  d=  7                       = 0.116
  d=  8                       = 0.530
  d= 12                       = 0.144
  d= 16                       = 0.530
```

Ratios match the Capelli prediction (powers of 2 lift the constraint
weight). No useful asymptotic compression visible.

### Q5 — exhaustive composite catch on `n ≤ 500`, `d ∈ {3,4,6}`

For `d=4, a=1`: every composite in `[4, 500]` caught at `b=1`. For
`d=3, a=1`: 4 composites need `b=2 or 3` (specifically `121` and the
prime cubes `9, 27, 81, 243`). No composite survives `b ∈ [1..11]`.

### Q6 — adversarial probe: 30 Carmichaels (561..410041), 13 base-2 PSPs, 11 prime squares

This is the empirically interesting result.

```
  --- d=3, a=1  (cyclotomic-suspect, x^3+1 = (x+1)*Phi_6) ---
    Carmichaels      19/30 pass at b=1     16/30 still pass for all b ∈ {1..4}
    Fermat-PSP base2  2/13                  2/13
    Prime squares    1/11                   0/11

  --- d=3, a=2  (Eisenstein, x^3+2 irreducible over Q) ---
    Carmichaels      0/30                   0/30
    Fermat-PSP base2 0/13                   0/13
    Prime squares    0/11                   0/11

  --- d=4, a=1  (x^4+1 = Phi_8, cyclotomic) ---
    Carmichaels      3/30 pass at b=1      2/30 still pass for b ∈ {1..4}
    PSP / squares     0                      0

  --- d=4, a=2  (non-cyclotomic) ---
    Carmichaels      0/30                   0/30
    PSP / squares    0                      0

  --- d=5, a=1 ---  Carm 2/30 pass        --- d=5, a=2 ---  0/30
  --- d=6, a=1 ---  Carm 5/30 pass        --- d=6, a=2 ---  0/30
  --- d=8, a=1 ---  Carm 0/30             --- d=8, a=2 ---  0/30
  --- d=12, a=1 --- Carm 0/30             --- d=12, a=2 --- 0/30
```

**Pattern:** at small `d`, the cyclotomic-flavour case `a = 1` admits a
non-trivial Carmichael leak (up to **19/30 at d=3**); the Eisenstein
case `a = 2` does NOT (`0/30` at every `d ∈ {3,4,5,6,8,12}` tested up
to `n = 410041`).

The Carmichaels that leak in the cyclotomic case are biased: e.g.
`1729 = 7·13·19`, `2821 = 7·13·31`, `8911 = 7·19·67`, `15841 =
7·31·73` all share small prime factor `7` and pass `d=3 a=1 b=1` —
their multiplicative orders mod `7^k` interact with the small
cyclotomic factor in a way `a=2` breaks.

### Q7 — scaling probe to `n ≤ 5000` (300 stratified samples)

```
  d=4, a=1, b=1:  primes 300/300, composites 0/300
  d=4, a=2, b=1:  primes 300/300, composites 0/300
  d=6, a=1, b=1:  primes 300/300, composites 0/300
```

Among 300 stratified composite samples up to 5000, no false positive
for any of the three configurations. The Carmichael leak in `a=1`
appears at the specific Carmichaels Q6 caught, not in the bulk.

## Interpretation

### Why `a=1` leaks Carmichaels and `a=2` does not

For `a = 1` and `d` admitting a small cyclotomic factor (`d ∈
{3,4,6,...}`), `x^d + 1` factors mod `n` over `(x+1)·Phi_{2d}(x)` or
similar. Korselt's criterion (Carmichael ⟺ `n` squarefree, `p-1 | n-1`
for every `p | n`) gives Carmichaels a *uniform* Frobenius-like
property mod every prime factor — and that property aligns with
small cyclotomic mod structure. So Carmichaels accidentally satisfy
the AKS-style identity in `(Z_n / (x^d + 1))*` even though they aren't
prime.

Setting `a = 2` (Eisenstein-flavour) destroys this alignment: `x^d +
2` has no rational roots and no cyclotomic structure, so the identity
becomes a "stronger" Frobenius-like check and the Carmichael
near-miss fails.

### Why this does NOT solve PRIMES in TC^0

Even with `a = 2` and small `d`, deciding the congruence

```
   (x + b)^n  =?  x^n + b   in   Z_n[x] / (x^d + a)
```

requires repeated squaring of an element of a `d`-dimensional `Z_n`
module — i.e. **growing-dim matrix powering** if we encode the squaring
as a `d × d` matrix multiplication. This is the very primitive E5.3
calls "the only OPEN frontier"; it is open at the TC^0 / NC^1
boundary. The choice of cyclotomic vs Eisenstein modulus changes the
*structure* of the matrix but not its *dimension*. So the depth
question is unchanged.

Furthermore, for the test to be a *theorem*-grade primality test (not
just an empirical filter), one needs an AKS-style correctness proof.
The standard AKS proof uses the cyclotomic structure exactly to embed
the introspective monoid into `F_p[x]/h(x)` of large enough degree.
Without that structure, the empirical absence of Carmichael leaks at
`n ≤ 410041` is encouraging but does **not** rule out larger
pseudoprimes — and there is no obvious analog of the AKS counting
argument for Eisenstein moduli.

### Net effect on Chain E

* No change to the open status of E5.3 (growing-dim MPOW in TC^0).
* Sub-attack 2 of FOCUS-1 is **closed**: a non-cyclotomic ring
  decomposition does not avoid the polylog-dim MPOW dependency.
* The empirical `a=1` vs `a=2` Carmichael separation is a small
  structural sharpening worth recording in CLOSED_PATHS but not
  worth its own `novel/` entry (it refines a known phenomenon — small
  cyclotomic leaks correlate with Korselt structure — rather than
  unveiling a new one).

## Reproducibility

All results from running this script with `python3
non_cyclotomic_ring.py` (≈ 60 seconds with sympy, no other deps).

## Next sub-attacks

* Sub-attack 1 (Bernstein 2003 smaller-r): independent line, still
  tractable.
* Sub-attack 3 (Healy-Viola Frobenius transplant): requires partial
  Frobenius analysis on `(Z/qZ)[x]/(x^r - 1)` for `q | n - 1` —
  independent of this construction, also still tractable.

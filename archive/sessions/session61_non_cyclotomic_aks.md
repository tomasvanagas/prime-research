# Session 61 — FOCUS-1 sub-attack 2: non-cyclotomic ring AKS

**Date:** 2026-04-26
**Mode:** normal (FOCUS-1, construction)
**Run state:** 36 → 37
**TODO item:** FOCUS-1 sub-attack 2 (non-cyclotomic ring decomposition)

## Summary

Built and probed an AKS-style primality congruence

```
   (x + b)^n  ≡  x^n + b    (mod  x^d + a, n)
```

with `f(x) = x^d + a` replacing the standard cyclotomic modulus
`x^r - 1`. Empirical sweep on `n` up to 410041 (incl. 30 Carmichaels
and 13 Fermat pseudoprimes base 2) and `d ∈ {3,…,12}`, `a ∈ {1,2}`.

## Headline result

**Closes FOCUS-1 sub-attack 2 as FAIL (E + C).** No bypass of E5.3
(growing-dim MPOW in TC^0). The Capelli irreducibility test is in TC^0
*conditional* on a prime factor of `n` being given — circular for
primality testing.

## Empirical observation

Clean separation between cyclotomic-flavour and Eisenstein modulus on
Carmichael adversarial cases:

| `d` | `x^d + 1` cyclotomic-flavour Carm leak | `x^d + 2` Eisenstein Carm leak |
|-----|----------------------------------------|--------------------------------|
| 3   | 19 / 30                                 | 0 / 30                         |
| 4   | 3 / 30                                  | 0 / 30                         |
| 5   | 2 / 30                                  | 0 / 30                         |
| 6   | 5 / 30                                  | 0 / 30                         |
| 8   | 0 / 30                                  | 0 / 30                         |
| 12  | 0 / 30                                  | 0 / 30                         |

Carmichaels sharing prime factor 7 (e.g. 1729=7·13·19, 2821=7·13·31,
8911=7·19·67, 15841=7·31·73) leak through small cyclotomic moduli of
the form `x^d + 1`. Korselt's `p - 1 | n - 1` for every `p | n`
interacts coherently with the small cyclotomic structure; the
Eisenstein twist `+a` for `a ≥ 2` breaks the alignment.

This is a small structural sharpening of CLOSED_PATHS line 690 (S47
cyclotomic-CRT closure) but is **not novel/-grade** — the underlying
Korselt/cyclotomic interaction is folklore in primality-testing
literature. No new `novel/` document.

## Why this does NOT bypass E5.3

* The congruence requires repeated squaring of an element of a
  `d`-dimensional `Z_n` module — that is, `d × d` matrix powering
  with `d = polylog(n)`. **Cyclotomic vs Eisenstein modulus changes
  the structure of the matrix but not its dimension.** The Mereghetti-
  Palano fixed-`k` MPOW result does not extend; the open frontier
  `growing-dim MPOW in TC^0?` (E5.3) is unaffected.
* No analog of the AKS introspective-monoid counting argument has
  been produced for Eisenstein moduli, so even if the `a = 2`
  Carmichael leak rate stays at 0 indefinitely, this is empirical
  evidence not theorem.
* Capelli's irreducibility test for `x^d + a` over `F_p` is itself
  in TC^0, but requires a prime factor `p | n` as input — circular
  for primality testing.

## Closures filed

* `status/CLOSED_PATHS.md` line 714 added: "Non-cyclotomic ring AKS
  variant Z_n[x]/(x^d+a) (FOCUS-1 sub-attack 2)" → FAIL, E+C.
* `status/SESSION_INSIGHTS.md` Session 61 entry added.
* No change to `OPEN_PROBLEMS.md` (Chain-E open frontier unchanged).
* No change to `novel/`.

## Files produced

* `experiments/circuit_complexity/aks_alternative/non_cyclotomic_ring/non_cyclotomic_ring.py`
* `experiments/circuit_complexity/aks_alternative/non_cyclotomic_ring/non_cyclotomic_ring_results.md`
* `archive/sessions/session61_non_cyclotomic_aks.md` (this file)

## Next steps for FOCUS-1

* Sub-attack 1 (Bernstein 2003 smaller-r AKS): independent attack
  surface, still untested. Read cr.yp.to/papers/aks.pdf for the
  strengthened gcd condition.
* Sub-attack 3 (Healy-Viola Frobenius transplant): independent
  attack surface, still untested. Compute power-Frobenius orbit
  structure on `(Z/qZ)[x]/(x^r - 1)` for small `(n, q, r)` triples
  with `q | n - 1`.

If both sub-attacks 1 and 3 close, Chain E becomes "computationally
cornered" and only Brandt MKtP (FOCUS-3) and a fundamentally new
lower-bound technique remain.

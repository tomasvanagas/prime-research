# Session 66 — FOCUS-1 sub-attack 1: Bernstein 2003 smaller-r AKS variant

**Date:** 2026-04-26
**Mode:** deep-focus task #1 (construction)
**Run state:** 43 → 44
**TODO item:** FOCUS-1 sub-attack 1 (Bernstein 2003 strengthened-gcd AKS)

## Summary

Built and probed the Bernstein-style strengthening of the AKS test:
the standard polynomial congruence

```
   (X + a)^n  ≡  X^n + a   (mod  n,  X^r - 1)
```

augmented with a deterministic gcd certificate on the residual
coefficients. The aim was to determine whether the strengthened
gcd condition (which Bernstein 2003 uses to deterministically
certify primality at the empirical small-r regime r ~ log²n)
admits a TC⁰-amenable formulation that would push PRIMES (or our
target π(x)) below the Chain-E open frontier `E5.3`.

## Headline result

**Closes FOCUS-1 sub-attack 1 with mode (E) Equivalence.**

The strengthened-gcd test is empirically very powerful — it
extracts a non-trivial factor of n in 13/13 composite failures
including all 7 sampled Carmichaels — but its complexity sits at
exactly the same NC¹/TC⁰ frontier as the growing-dim MPOW
primitive `E5.3` it purports to replace. With this closure, all
three FOCUS-1 sub-attacks are now closed and **Chain E is
computationally cornered** per `TODO.md`.

## Concrete numbers (S47-aligned grid, 22 samples in [101, 10⁶])

| Quantity | Result |
|----------|--------|
| Mean `r / log²(n)` | 1.207 |
| Mean `r − log²(n)` | 25.05 |
| `r` is prime | 21/22 |
| Mean `max_dim/r` (cyclotomic CRT) | 0.9899 (= S47) |
| Primes pass canonical-r test | 3/3 |
| Composites pass canonical-r test | 0/13 |
| Composites yielding gcd-extracted factor | **13/13** |
| Carmichaels caught with factor | 7/7 |

## Why the gcd extraction succeeds empirically

For prime n, `(X+a)^n = X^n + a^n = X^n + a` in `Z_n[X]/(X^r − 1)`
because (i) Frobenius `(X+a)^n ≡ X^n + a^n mod n` and (ii) Fermat
`a^n ≡ a mod n`. For composite n, both fail:

* The middle binomials `binom(n, k) mod n` are non-zero for k a
  proper divisor of any prime factor of n,
* These appear as coefficients of the residual `(X+a)^n − (X^n + a)`,
* Each residual coefficient is `≡ a^{n−k} · binom(n, k) mod n`,
* `gcd(binom(n,k), n)` is a non-trivial factor of n whenever k is
  a "small base-p sub-string" of n for some prime p | n
  (Lucas's theorem applied prime-by-prime via CRT on Z_n).

So composites leak factors through their binomial structure. This
is the *same* mechanism as S64's residual analysis (mod-q Frobenius
transplant) but operates over Z_n directly — preserving the AKS
Frobenius identity primes do satisfy, rather than reducing it away.

## Why the test, even working, would not be TC⁰

The strengthened-gcd test requires:

* Length-`O(log n)` sequential repeated squaring of an r-dim
  polynomial in Z_n[X]/(X^r − 1). Each squaring is an `r × r`
  linear map over Z_n. The chain is `O(log n)` deep.
* For each of `S = √φ(r) · log n` values of a, compute up to `r`
  integer gcds of `O(log n)`-bit coefficients with n.

| Primitive | NC¹ | TC⁰ |
|---|---|---|
| growing-dim r × r MPOW over Z_n | YES | OPEN (E5.3) |
| gcd of two ℓ-bit integers | YES (Hesse-Allender-Barrington 2002) | OPEN, conjectured NO |
| polynomial multiplication mod (n, X^r − 1) | YES | OPEN |

The strengthening **trades** one NC¹/TC⁰-frontier problem for
another. Both pieces are in NC¹. Neither is known in TC⁰. Even
if integer gcd were resolved in TC⁰ tomorrow, the polynomial-
multiplication step would still be growing-dim r × r MPOW over
Z_n — exactly E5.3.

This is a clean (E) closure — the depth question is preserved.

## Why we don't elevate this to a `novel/` entry

The factor-extraction observation (gcd of residual coefficient
yields a divisor of composite n in 13/13 cases including all 7
Carmichaels) is empirically striking but follows from the
standard AKS proof structure: residual coefficients carry
multiplicative-character information from binom(n, k) which
carries n's prime structure by Lucas+CRT. Both AKS variants
(Bernstein 2003, Lenstra-Pomerance 2005) and the surrounding
literature implicitly use this. It is not in `novel/` because:

1. It does not unblock any complexity question.
2. It is implicit in the standard AKS proof.
3. It is at most an algorithmic curiosity (AKS-as-side-channel-
   factoring, an old folklore observation).

The S66 contribution is the explicit empirical verification on
the S47 grid plus the precise NC¹/TC⁰ frontier-equivalence
argument that closes sub-attack 1.

## Closures filed

* `status/CLOSED_PATHS.md` line ~722 added: "Bernstein 2003
  strengthened-gcd AKS variant (FOCUS-1 sub-attack 1)" → FAIL,
  E (with weak C).
* `status/OPEN_PROBLEMS.md` updated with the S66 status note —
  Chain E is computationally cornered.
* `status/SESSION_INSIGHTS.md` — Session 66 entry added.
* `TODO.md` FOCUS-1 status updated to reflect all three sub-
  attacks closed.

## What's left after S66

The AKS family of TC⁰ approaches is exhausted. The remaining
levers on the only open direction (Problem #1 in
`OPEN_PROBLEMS.md`) are:

1. **Brandt MKtP** (FOCUS-3) — un-engaged framework, possible
   diagonalization argument for natural functions; `EDGES.md`
   notes no follow-up papers since 2024 TCC.
2. **Fundamentally new lower-bound technique** — open problem
   class, no concrete proposal.
3. **Non-AKS TC⁰ primality test** — using only scalar operations,
   long-standing aspiration since S15. No active line of attack.
4. **Engineering improvements** to `algorithms/v10_c_accelerated.py`
   (Gourdon variant, segmented sieve, SIMD).
5. **Literature watch** (FOCUS-5).

## Files produced

* `experiments/circuit_complexity/aks_alternative/bernstein_smaller_r/bernstein_smaller_r.py` (~9 KB)
* `experiments/circuit_complexity/aks_alternative/bernstein_smaller_r/bernstein_smaller_r_results.md` (~7 KB)
* `archive/sessions/session66_bernstein_smaller_r.md` (this file)

## Cross-references

* `EDGES.md` E5.3 — open Chain-E frontier (growing-dim MPOW in TC⁰?)
* `status/CLOSED_PATHS.md` lines 690 (S47), 714 (S61), 719 (S64),
  ~722 (S66) — the four FOCUS-1 closures
* `TODO.md` FOCUS-1 — all three sub-attacks closed

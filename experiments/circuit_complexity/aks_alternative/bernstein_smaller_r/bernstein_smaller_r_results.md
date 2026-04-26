# Bernstein 2003 strengthened-gcd AKS variant — results

**Construction (S66 / Run 43, FOCUS-1 sub-attack 1 from `TODO.md`).**
Tests whether the smaller-r AKS variant due to Bernstein 2003 (`Proving
primality after Agrawal-Kayal-Saxena`, cr.yp.to/papers/aks.pdf) — in
which the standard order-condition `ord_r(n) > log²(n)` is augmented
with a stronger gcd condition that *deterministically* certifies
primality at smaller r — admits a TC⁰-amenable formulation. With
sub-attacks 2 (S61, line 714) and 3 (S64, line 719) already closed,
this is the only remaining un-built piece of FOCUS-1.

## Verdict

**FAIL — closes FOCUS-1 sub-attack 1 with closure mode (E) Equivalence.**

The Bernstein-style strengthening *does* discriminate composites with
striking empirical power on the S47 grid (gcd extracts a non-trivial
factor of n in **13/13 composite failures, including all 7 Carmichaels
tested**), but the strengthened gcd condition itself is in **NC¹** —
the same NC¹/TC⁰ frontier as the growing-dim MPOW primitive `E5.3` it
purports to replace. No depth collapse. No reduction in the matrix-
powering dimension `r-1` either.

With this closure, **all three FOCUS-1 sub-attacks are now closed**;
Chain E should be marked "computationally cornered" per `TODO.md`.

## Q1 — Empirical r-values on S47 22-sample n grid

```
         n |      r |  log2(n)^2 |  r-log^2 |  max_dim |  r prime
       101 |     41 |         36 |        5 |       40 |     True
       561 |     89 |         81 |        8 |       88 |     True
      1009 |     83 |         81 |        2 |       82 |     True
      1024 |    227 |        100 |      127 |      226 |     True
      1729 |    107 |        100 |        7 |      106 |     True
      4096 |    311 |        144 |      167 |      310 |     True
     10007 |    173 |        169 |        4 |      172 |     True
     65537 |    271 |        256 |       15 |      270 |     True
    100003 |    263 |        256 |        7 |      262 |     True
    524287 |    331 |        324 |        7 |      330 |     True
   1000003 |    379 |        361 |       18 |      378 |     True
        ...  (22 rows total — full table in run log)
```

**Headline:** mean `r/log²(n) = 1.207`, mean `(r − log²n) = 25.05`,
`r` is prime in **21/22** cases.

The empirical `r` from the standard AKS order-condition is already
within an additive-constant of the Bernstein conjectured bound
`r = O(log² n)`. *"Smaller r"* is therefore not the gap Bernstein
2003 closes — small `r` is empirical reality already. What Bernstein
2003 attempts to provide is a **deterministic correctness theorem**
at this r, via the strengthened gcd condition.

The two outliers (`r/log²n` ≈ 2.27 for n=1024 and 2.16 for n=4096)
are powers of 2: every odd `r ∈ {3,..,225}` happens to satisfy
`gcd(r, 2^k) = 1` but the order condition needs the multiplicative
order of 2 mod `r` to be large, and 2 has small order modulo many
small primes (e.g. ord₅(2) = 4, ord₇(2) = 3). A non-issue: powers of
2 are trivially composite by inspection; no algorithmic significance.

## Q2/Q3 — AKS polynomial congruence and Bernstein-gcd extraction

```
         n |     r | is_prime(n) |  pass | fails |     gcd_witness
       101 |    41 |        True |  True |     0 |            None
       102 |    59 |       False | False |     1 |     (1, 51, 51)
       561 |    89 |       False | False |     1 |   (1, 374, 187)
      1001 |    83 |       False | False |     1 |     (1, 777, 7)
      1009 |    83 |        True |  True |     0 |            None
      1024 |   227 |       False | False |     1 |   (1, 256, 256)
      1105 |   103 |       False | False |     1 |   (1, 884, 221)
      1729 |   107 |       False | False |     1 |   (1, 266, 133)
      2048 |   131 |       False | False |     1 |   (1, 768, 256)
      2465 |   139 |       False | False |     1 |   (1, 1189, 29)
      2821 |   149 |       False | False |     1 |  (1, 1085, 217)
      4096 |   311 |       False | False |     1 |  (1, 1536, 512)
      6601 |   151 |       False | False |     1 |  (1, 3059, 161)
      8911 |   179 |       False | False |     1 | (1, 5092, 1273)
     10001 |   197 |       False | False |     1 |   (1, 548, 137)
     10007 |   173 |        True |  True |     0 |            None
```

`gcd_witness = (a, c, g)` means: with the smallest discriminating `a`,
the residual `(X+a)^n − (X^n + a) mod (n, X^r−1)` has a non-zero
coefficient `c` such that `g = gcd(c, n)` is a non-trivial factor of n.

**Discrimination:** primes 3/3 pass, composites 0/13 pass.

**Bernstein-gcd extraction:** **13/13 composites yield a non-trivial
factor of n**. Every Carmichael in the sample (561, 1105, 1729, 2465,
2821, 6601, 8911) gives a factor:

| n | factorisation | gcd witness | extracted factor |
|---|---|---|---|
| 561 | 3·11·17 | (1, 374, 187) | 11·17 = 187 |
| 1105 | 5·13·17 | (1, 884, 221) | 13·17 = 221 |
| 1729 | 7·13·19 | (1, 266, 133) | 7·19 = 133 |
| 2465 | 5·17·29 | (1, 1189, 29) | 29 |
| 2821 | 7·13·31 | (1, 1085, 217) | 7·31 = 217 |
| 6601 | 7·23·41 | (1, 3059, 161) | 7·23 = 161 |
| 8911 | 7·19·67 | (1, 5092, 1273) | 19·67 = 1273 |

This is empirically much stronger than I expected before running. The
residual coefficient `c` is congruent to a polynomial in `binom(n, k)`
values for various `k`, and Lucas's theorem gives these binomials a
multiplicative-character structure mod each prime divisor of n that
collides predictably for Carmichael-like `n`. The extraction is
identical in spirit to the residual analysis in S64 (Frobenius
transplant) but operates over ℤ_n directly rather than via a mod-q
twist — and unlike S64 which lost information by reducing mod q, this
*preserves* the Frobenius identity that primes do satisfy.

**Mild novelty note (Q3):** the AKS-with-gcd-extraction algorithm
simultaneously certifies primality AND, on composite input, returns
a non-trivial factor of n. I have not seen this stated explicitly in
the literature on AKS variants. It is a corollary of the standard
AKS argument — the residual of (X+a)^n − (X^n + a) is a polynomial
whose coefficients are congruent to `a^{n−k}·binom(n, k)` and these
binomials are non-zero mod n iff `n` is composite — but it is worth
recording as an observation. It does not constitute a `novel/` entry
because (a) it does not unblock any complexity question, and (b) it
is implicit in any AKS proof that explains why the test works.

## Q4 — dim/r ratio (Bernstein test inherits the same matrix-powering dim)

```
         n |     r | max cyc dim |  dim/r | r prime?
       101 |    41 |          40 |  0.976 |     True
      1024 |   227 |         226 |  0.996 |     True
      4096 |   311 |         310 |  0.997 |     True
    100001 |   289 |         272 |  0.941 |    False    <- only composite r
   1000003 |   379 |         378 |  0.997 |     True
```

**Mean `max_dim/r` = 0.9899** across all 22 samples — identical to
S47. The CRT decomposition of `x^r − 1` saves at most a factor of 1.06
on dimension, and only when `r` happens to be composite. Bernstein's
strengthening operates on the *gcd-test side*, not the matrix
dimension; it inherits whatever dim the underlying ring has.

## Q5 — small-r probe (do composites leak at r below the order bound?)

```
       n | is_prime | r_AKS | log_n+1 |               r_test results
     101 |     True |    41 |       7 |    r=7:P, r=12:P, r=17:P, r=20:P, r=40:P
     102 |    False |    59 |       7 |                            r=7:F, r=29:F
     561 |    False |    89 |      10 |                                   r=10:F
    1009 |     True |    83 |      10 |   r=10:P, r=18:P, r=31:P, r=41:P, r=82:P
    1024 |    False |   227 |      11 |                  r=11:F, r=31:F, r=113:F
```

Surprise: even at `r ≪ r_AKS` (e.g. `r = log₂(n) + 1`), the polynomial
test still discriminates the *specific* composites in our sample. This
is because (a) prime-power composites like `1024` violate the test
trivially via `binom(2^10, k)` divisibility patterns, and (b) Carmichael
561 has small odd-prime factors that already break the cyclotomic
identity at small r. The AKS bound `r ~ log² n` is a *worst-case* bound
needed to defeat *adversarially constructed* near-primes; typical
composites are caught at far smaller r.

This is an empirical complementary observation to Q3: even pre-Bernstein
AKS would work at `r = O(log n)` for *most* composites, but cannot be
proven correct without either the order-condition bound or some
auxiliary check. Bernstein 2003 supplies the latter.

## Q6 — TC⁰ placement of the strengthened gcd condition

The strengthened condition that lets AKS run at empirical r ~ log²(n)
is a gcd test on the residual polynomial coefficients:

> For each `a ∈ [1, S]`, compute `d_a(X) = (X+a)^n − (X^n + a)`
> in `ℤ_n[X]/(X^r − 1)`. For every non-zero coefficient `c` of `d_a`,
> check `gcd(c, n)`. Test passes iff every `d_a` is identically 0.

The polynomial multiplications produce coefficients in ℤ_n, i.e. each
O(log n) bits. The gcd is on two integers of O(log n) bits each.

| Primitive | NC¹? | TC⁰? | Source |
|---|---|---|---|
| gcd of ℓ-bit integers | YES | OPEN, conjectured NO | Hesse-Allender-Barrington 2002 |
| polynomial mult mod (n, X^r−1) | YES | OPEN | E5.3 (growing-dim MPOW) |
| binary integer mod | YES | YES | Hesse 2001 |

The strengthening **replaces one NC¹/TC⁰-frontier problem (growing-dim
MPOW) with another (gcd of polylog-bit integers, repeated S = √φ(r) ·
log n times)**. Both problems sit at the same depth frontier. No
progress on placing PRIMES in TC⁰ results from the substitution.

Crucially: even if integer gcd were placed in TC⁰ tomorrow (a major
result), the polynomial-multiplication step in Bernstein's test is
itself a growing-dim r×r linear map over ℤ_n — exactly E5.3. So the
gcd strengthening is at best a re-arrangement of the same primitive,
at worst a strict equivalence; never an escape.

## Q7 — Why this closes Chain E sub-attack 1

| Failure mode | Argument |
|--------------|----------|
| **(E) Equivalence** | The strengthened-gcd test does length-O(log n) growing-dim r×r MPOW over ℤ_n exactly as standard AKS, then performs S = √φ(r)·log n integer-gcd checks of O(log n)-bit integers. Both pieces are NC¹-known / TC⁰-open. The substitution has NO effect on the depth question. |
| **(C) Circularity (weak)** | If we abandon the order-condition r-bound and try to set `r = log²(n)` deterministically, we need a sufficient certificate that the test is correct at this r. Bernstein's certificate is the gcd extraction. Verifying that no Carmichael-like adversary can simultaneously satisfy the polynomial congruence and have all residual coefficients gcd-trivially-coprime to n is itself an open problem requiring AKS-style introspective-monoid counting. The certificate doesn't avoid the proof obligation; it relocates it. |

## Closures filed

* `status/CLOSED_PATHS.md`: append "Bernstein 2003 strengthened-gcd
  AKS variant (FOCUS-1 sub-attack 1)" → FAIL, E (with weak C),
  session 66.
* `status/SESSION_INSIGHTS.md`: append session 66 entry.
* `status/OPEN_PROBLEMS.md`: mark Chain E as **computationally
  cornered** per `TODO.md` FOCUS-1 directive (all three sub-attacks
  closed, only Brandt MKtP / new lower-bound technique remain).
* `TODO.md`: update FOCUS-1 status: sub-attack 1 closed, mark
  computationally-cornered milestone reached.

## State of FOCUS-1 after Session 66

| Sub-attack | Status |
|-----------|--------|
| 1 — Bernstein 2003 strengthened gcd | **CLOSED** (S66, this) — (E) Equivalence |
| 2 — Non-cyclotomic ring AKS Z_n[x]/(x^d+a) | CLOSED (S61, line 714) — (E)+(C) |
| 3 — Healy-Viola Frobenius transplant | CLOSED (S64, line 719) — (E)+(I) |

**All three sub-attacks closed.** Chain E (TC⁰ via growing-dim
matrix powering) is now "computationally cornered": every modulus-
twist and gcd-strengthening of the AKS family reduces to growing-dim
MPOW over the same r ~ log² n. The only remaining escapes are:

1. **Brandt MKtP** (FOCUS-3 in TODO.md) — un-engaged framework with
   a possible diagonalization argument for natural functions.
2. **A fundamentally new lower-bound technique** (open problem class).
3. **A non-AKS TC⁰ primality test** using only scalar operations
   (mentioned as "what might still work" in OPEN_PROBLEMS.md).

The first two are exploratory and outside the AKS family; the third
has been an unresolved aspiration since S15.

## Files

* `bernstein_smaller_r.py` (≈8 KB) — construction and tests
* `bernstein_smaller_r_results.md` — this file

## Cross-references

* `EDGES.md` `E5.3` — open Chain-E frontier (growing-dim MPOW in TC⁰?)
* `status/CLOSED_PATHS.md` line 690 (S47, cyclotomic-CRT), line 714
  (S61, non-cyclotomic ring), line 719 (S64, Frobenius transplant) —
  adjacent FOCUS-1 closures
* `TODO.md` FOCUS-1 sub-attack 1 — closed by this experiment
* `archive/sessions/session66_bernstein_smaller_r.md` — session
  synthesis

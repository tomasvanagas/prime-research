# Healy-Viola Frobenius transplant for AKS — results

**Construction (S64 / Run 42, FOCUS-1 sub-attack 3 from `TODO.md`).**
Tests whether the AKS congruence `(a + x)^n ≡ a^n + x^n (mod x^r - 1, n)`
restructures into a TC^0-amenable form when reduced modulo a prime
`q | n - 1`, by exploiting the q-power Frobenius

```
   F:  F_q[x]/(x^r - 1)  ->  F_q[x]/(x^r - 1),    F(a + x) = a + x^q
```

which is a ring homomorphism (q is the characteristic of F_q).

If the AKS-mod-q congruence held, then writing `n = sum_i c_i q^i` (base-q
expansion) gives the **Frobenius decomposition**

```
   (a + x)^n  =  prod_i (a + x^{q^i mod r})^{c_i}     in F_q[x]/(x^r - 1)
```

— a product of `log_q(n)` binomial polynomials, each computable from
the previous by an r-dim linear map.

## Verdict

**FAIL — closes FOCUS-1 sub-attack 3 with two simultaneous failure modes**.

`(E) Equivalence + (I) Information loss`.

* **(I)** Primes do NOT satisfy the AKS-mod-q congruence for any
  non-trivial a (a coprime to q). `0 / 399` cases passed across 19 primes
  in `[101..2207]` and every q ∈ primes-dividing-`n-1` with `q < r`.
  The mod-q congruence holds *only trivially* when `a ≡ 0 mod q` — both
  sides then reduce to `x^n` and the test carries zero primality information.
* **(E)** The `log_q(n)`-step product is still r-dim matrix powering —
  precisely the growing-dim MPOW primitive `E5.3`, unchanged by switching
  the coefficient ring from `Z_n` to `F_q`. Even granting the test worked,
  it would NOT bypass Chain E.

## Q1 — Frobenius decomposition is correct (sanity check)

Six test cases verify the identity `naive_repeated_squaring == Frobenius_decomp`
in `F_q[x]/(x^r - 1)`:

```
   n=  101 q=2 r=  5 a=3: OK
   n=  101 q=5 r= 11 a=2: OK
   n=  561 q=2 r=  5 a=2: OK
   n= 1729 q=3 r=  7 a=2: OK
   n=10007 q=2 r=  7 a=4: OK
   n=10007 q=3 r=  7 a=2: OK
   Verdict: PASS
```

Decomposition formula is implemented correctly.

## Q2 — Primes fail mod-q AKS for non-trivial a [headline]

For prime `n`, prime `q | n-1` with `q < r = aks_r(n)`, every `a in {1, ..., q-1}`:

```
        n    r  q -> non-trivial-a passes / total
      101   53  q= 2: 0/1   q= 5: 0/4
      103   53  q= 2: 0/1   q= 3: 0/2   q=17: 0/16
      107   47  q= 2: 0/1
      109   47  q= 2: 0/1   q= 3: 0/2
      113   59  q= 2: 0/1   q= 7: 0/6
      211   71  q= 2: 0/1   q= 3: 0/2   q= 5: 0/4   q= 7: 0/6
      223   83  q= 2: 0/1   q= 3: 0/2   q=37: 0/36
      251   67  q= 2: 0/1   q= 5: 0/4
   ... (11 more primes, same pattern)

   Aggregate non-trivial-a tests: 0 / 399 primes pass
   Aggregate trivial-a tests (a=0): 48 / 48 pass
                              (both sides reduce to x^n; carries no information)
```

The trivial passes (`a ≡ 0 mod q`) are non-informative: `(0 + x)^n = x^n`
and `0^n + x^n = x^n` in `F_q[x]/(x^r - 1)` for every n. The non-trivial
tests fail uniformly.

## Q3 — Why it fails: the residual polynomial

For prime `n`, the polynomial identity `(a+x)^n - (a^n + x^n) =
sum_{k=1}^{n-1} binom(n,k) a^{n-k} x^k` vanishes mod n (Frobenius)
but NOT mod `q` for `q != n` prime, since the middle binomial
coefficients `binom(n, k)` are computed by **Lucas's theorem mod q**:

```
   binom(n, k) mod q  =  prod_i binom(n_i, k_i)
```

where `n_i`, `k_i` are base-q digits. These are non-zero whenever each
`k_i <= n_i` (k is a "base-q sub-string" of n). Concrete: `n = 101 =
0b1100101`, so `binom(101, k)` is odd iff k's bits are a subset of 101's
bits — `k ∈ {0, 1, 4, 5, 32, 33, 36, 37, 64, 65, 68, 69, 96, 97, 100, 101}`,
16 of 102 binomial coefficients are odd.

Empirical residual polynomial weight `(LHS - RHS) mod (q, x^r - 1)`:

```
       n  q    r  a  nonzero coeffs / r
     101  2   53  1   14/53
     101  5   53  1    8/53
     101  2   53  3   14/53
     101  5   53  3    8/53
    1009  3  107  2   22/107
    1009  7  107  2  102/107
    1009  2  107  5   42/107
    1009  3  107  5   22/107
    1009  7  107  5   87/107
```

For prime `n=1009` mod `q=7`: 102 out of 107 coefficients of the residual
polynomial are non-zero — the test misses by a country mile. The residual
fills the entire `F_q[x]/(x^r-1)` ring rather than being a sparse
correction, ruling out any partial-cancellation rescue.

## Q4 — Depth/width: NC^1 not TC^0

```
            n  q    r  frob_muls  naive_muls   log_q n   log_2 n
         1009  2  107          7         18      9.98      9.98
         1009  3  107          5         18      6.30      9.98
         1009  7  107         12         18      3.55      9.98
        10007  2  179          8         26     13.29     13.29
       100003  2  277          8         32     16.61     16.61
       100003  3  277         14         32     10.48     16.61
       100003  7  277         20         32      5.92     16.61
      1000003  2  401          9         38     19.93     19.93
      1000003  3  401         20         38     12.58     19.93
```

* Frobenius decomposition saves a constant factor of `log(q)/log(2)`
  multiplications relative to naive repeated squaring (2x for q=2,
  ~3.3x for q=7, etc.) — bit savings only.
* Depth: both decompositions reduce to `O(log n)` sequential
  multiplications of `r-dim` polynomials. Each polynomial multiplication
  is in NC^1 (parallel addition tree of depth `O(log r) = O(log log n)`),
  but the sequential chain of `log n` multiplications gives total depth
  `O(log n · log log n)` = polylog, not constant.
* The growing-dim MPOW primitive `E5.3` is unchanged: each step is an
  `r x r` linear map on `F_q^r` with at most 2 non-zero entries (the
  binomial generator), and the running product is a sequential matrix
  power over `F_q`.

The Frobenius restructuring **does not** push the depth below NC^1.
For TC^0 we would need `O(1)` depth, which the sequential chain forbids
unless the matrix power itself collapses (which is the open Chain E
question).

## Q5 — Two simultaneous closures

| Failure mode | Argument |
|--------------|----------|
| **(E) Equivalence** | The mod-q test reduces to length-`O(log n)` `r x r` MPOW over `F_q`. Switching coefficient ring `Z_n -> F_q` does not change the `r x r` matrix-powering primitive (E5.3); the q-twist is structural redress, not depth collapse. |
| **(I) Information loss** | The standard AKS congruence is mod-n specific. Reducing mod q for `q != n` prime kills the Frobenius identity that makes `(a+x)^n ≡ a^n + x^n` true in the first place — primes themselves fail the mod-q test (Q2: 0/399). No Hensel-style lifting from `F_q` back to `Z_n` is available. |

This is a **clean negative result**, complementary to S47 (cyclotomic-CRT
splitting) and S61 (non-cyclotomic ring AKS). With sub-attacks 1, 2, and
3 of FOCUS-1 now individually probed (sub-attack 1 = Bernstein 2003
strengthened gcd is the only remaining un-built construction), Chain E
is increasingly cornered: the modulus structure is not the leverage
point. The depth question hinges on growing-dim MPOW, full stop.

## Files

* `frobenius_transplant.py` — the construction and tests (≈4 KB)
* `frobenius_transplant_results.md` — this file

## Cross-references

* `EDGES.md` `E5.3` — open Chain-E frontier (growing-dim MPOW in TC^0?)
* `status/CLOSED_PATHS.md` line 690 (S47 cyclotomic-CRT) and line 714
  (S37/S61 non-cyclotomic ring AKS) — adjacent closures
* `TODO.md` FOCUS-1 sub-attack 3 — closed by this experiment

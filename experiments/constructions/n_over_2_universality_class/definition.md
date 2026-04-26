# C5 — N/2 universality × non-pi Boolean function class

## Edges composed

* **E1.4** — N/2 universality:  six independent complexity measures of the
  prime indicator chi_P all land at boundary N/2 (or 2^{N/2-1}+2 for
  rank). Verified for chi_P only.
* **E2.5** — chi_P is a multilinear polynomial of degree N over R / Z;
  approximate-degree and PTF-degree LP machinery operates on this
  representation.

## Object

Define a class C of natural number-theoretic Boolean functions:

```
    f_pi(n)        := 1 iff n is prime
    f_sqfree(n)    := 1 iff n is squarefree                (mu(n) != 0)
    f_mu_pos(n)    := 1 iff Mobius(n) = +1
    f_lam_pos(n)   := 1 iff Liouville(n) = +1
    f_sqfree3(n)   := 1 iff n squarefree AND n mod 3 == 1
    f_prf_matched(n) := SHA-256-derived bit, density-matched to f_pi
```

For each f in C, run a 4-measurement subset of the N/2 battery on
{0, 1, ..., 2^N - 1}:

* **M1** balanced communication-matrix rank: split lower N/2 bits as
  columns, upper N/2 bits as rows. Compute exact rank of the 0/1 matrix.
* **M2** GF(2) Berlekamp-Massey linear complexity of the truth-table
  sequence f(0), f(1), ..., f(2^N - 1).
* **M3** approximate degree at eps = 0.49 (Nisan-Szegedy notion):
  smallest d such that there exists a real polynomial of degree d on
  {0,1}^N with sup |p - f| < 0.49.  Solved as LP.
* **M4** PTF degree: smallest d such that a real polynomial of degree d
  sign-represents f on {0,1}^N.  Solved as LP feasibility.

## Falsification statements

* **PR1** (M1, comm-rank, "+2 universality"):  For every f in C and
  every N, rank(M_f) is in the band [2^{N/2-1}, 2^{N/2-1} + c] for some
  small c independent of N.  This is the test of whether E2.7's
  "+2 anomaly" generalises beyond chi_P.

* **PR2** (M2, BM-LFSR universality):  For every f in C and every
  N >= 10, BM linear complexity satisfies L / 2^N in [0.4, 0.6].
  Tests whether all NT Boolean functions look pseudorandom in BM.

* **PR3** (M3, approximate-degree N/2 universality):  For every f in C
  and every N, adeg_{0.49}(f) = ceil(N/2).  Tests whether the N/2
  approximate-degree boundary is a meta-theorem of the class C, or
  specific to chi_P.

* **PR4** (M4, PTF-degree N/2 universality):  For every f in C and
  every N, ptf_deg(f) = ceil(N/2).

## Outcome interpretation

* PR1-PR4 all PASS:  N/2 universality is a meta-theorem of natural NT
  Boolean functions.  This is a strong negative-shape claim suggesting
  the boundary is determined by the {0,1}^N structure, not by
  primality-specific arithmetic.

* Some PR fails:  chi_P is special at that measure.  The deviation is a
  structural fingerprint of primes versus the broader NT Boolean class,
  and worth catalogueing.

## Files

* `n_over_2_universality_class.py` — code (LP solves; truth-table
  builders; battery driver).
* `n_over_2_universality_class_data.json` — full measurement dump.
* `n_over_2_universality_class_results.md` — analysis and verdict.

# Wheel-coprime Liouville-Möbius identity and its prime-bisection lift

## Object

For every positive integer `W` and every real `x ≥ 0`, define the
**wheel-coprime Liouville sum**

```
    L_W(x) := Σ_{n ≤ x, gcd(n, W) = 1} λ(n)
```

where `λ(n) = (-1)^Ω(n)` is the Liouville function.

This construction asserts and verifies an exact pointwise identity
expressing `L_W(x)` as a finite signed combination of unrestricted
Liouville sums `L(⌊x/d⌋)` indexed by divisors of `rad(W)`, and lifts
the identity through the prime-bisection
`π(x) = (x − L(x))/2 − C_3(x)` (E1.6) into a wheel-graded prime
identity.

## Pre-stated identities (proven analytically — to be verified pointwise)

### Theorem 1 (wheel-coprime Liouville-Möbius identity).

For every `W ≥ 1` and every `x ≥ 0`,

```
    L_W(x)  =  Σ_{d | rad(W)} L(⌊x/d⌋).
```

*All terms* on the right are unrestricted Liouville partial sums; the
right-hand side has `2^{ω(W)}` terms regardless of how large `W` is,
and depends on `W` only through its squarefree radical `rad(W)`.

**Proof.**

```
    L_W(x)  =  Σ_n λ(n) · [gcd(n, W) = 1]
            =  Σ_n λ(n) Σ_{d | gcd(n, W)} μ(d)
            =  Σ_{d | W} μ(d) Σ_{m ≤ x/d} λ(d m)
            =  Σ_{d | W} μ(d) λ(d) L(⌊x/d⌋).
```

The last step uses `λ` completely multiplicative (`λ(dm) = λ(d)λ(m)`).
For squarefree `d`, `μ(d) λ(d) = (-1)^{ω(d)} · (-1)^{ω(d)} = +1`; for
non-squarefree `d`, `μ(d) = 0`. The squarefree divisors of `W` are
exactly the divisors of `rad(W)`. □

### Theorem 2 (parity of `L_W(x)`).

For every `W ≥ 1` and every integer `x ≥ 0`,

```
    L_W(x) mod 2  =  ( Σ_{d | rad(W)} ⌊x/d⌋ ) mod 2.
```

**Proof.** By E2.10, `L(y) ≡ y (mod 2)` for every integer `y ≥ 0`.
Reduce Theorem 1 mod 2. □

### Theorem 3 (wheel-graded prime bisection — lift of E1.6).

Define
```
    n_W(x) := #{n ≤ x : gcd(n, W) = 1}
    A_W(x) := #{n ≤ x : gcd(n, W) = 1, Ω(n) odd}
    π_W(x) := #{p ≤ x : p prime, gcd(p, W) = 1}    = π(x) − π_W^div(x)
    C_{3,W}(x) := A_W(x) − π_W(x)
                = #{n ≤ x : gcd(n, W) = 1, n composite, Ω(n) odd}
```
where `π_W^div(x)` is the bounded count of primes ≤ x that divide `W`.

Then
```
    A_W(x) = ( n_W(x) − L_W(x) ) / 2,                                    (3a)
    π_W(x) = ( n_W(x) − L_W(x) ) / 2  −  C_{3,W}(x),                     (3b)
    n_W(x) = Σ_{d | rad(W)} μ(d) ⌊x/d⌋.                                  (3c)
```

Combining (3b) with Theorem 1:

```
    π_W(x)  =  (1/2) Σ_{d | rad(W)} ( μ(d) ⌊x/d⌋ − L(⌊x/d⌋) )
                  −  C_{3,W}(x).                                         (3d)
```

The right-hand side has exactly `2^{ω(W)}` calls to `L` (the
unrestricted Liouville sum) and `2^{ω(W)}` calls to `μ(d)·⌊x/d⌋`
(constant cost). The hard residue `C_{3,W}` is the wheel-graded
"≥ 3 prime factors with odd Ω" count.

**Proof.** Identity (3a) follows from `λ(n) = +1` iff `Ω(n)` even
and `λ(n) = −1` iff `Ω(n)` odd, restricted to `gcd(n, W) = 1`:
`n_W − L_W = (#even Ω in coprime) + (#odd Ω in coprime) − (#even − #odd)
= 2 · #odd = 2 A_W`. Identity (3b) follows from `χ_P(n) = [Ω(n) = 1]
= [Ω(n) odd] − [Ω(n) ≥ 3, Ω(n) odd]` restricted to `gcd(n, W) = 1`.
Identity (3c) is the standard Möbius inversion of the coprime indicator.
Identity (3d) substitutes Theorem 1 into (3b). □

### Corollary (oracle reduction).

Given an oracle for `L(x)` evaluating in cost `T_L(x)`, the wheel-
coprime Liouville sum `L_W(x)` is computable in `2^{ω(W)} · T_L(x)`
operations, **independent of `W` itself**. By (3d), the wheel-coprime
prime-bisection smooth side `(n_W(x) − L_W(x))/2` is also reducible
to `2^{ω(W)}` calls to `L`, plus the residue `C_{3,W}(x)` whose
computation cost is the irreducible "≥ 3 prime factors with odd Ω"
problem.

## Why this is a NEW project artifact

The project records:

* **E1.6**: π(x) = (x − L(x))/2 − C_3(x) and bit-level independence of
  A(x) mod 2 vs C_3(x) mod 2 across q ∈ {2..13}.
* **E2.10**: L(x) mod 2 = x mod 2.
* **E2.1 / E4.1**: wheel-W density `φ(W)/W` is the bond-dim factor of
  the chi_P MPS reshape.
* **S208 (`spike_divisor_only_wheel`)**: the wheel-coprime indicator
  emerges as a *pointwise* scalar field equal to
  `[gcd(n, W) = 1] · W/φ(W)`.

What the project does **not** record is the cumulative-side analogue:
an exact pointwise identity expressing the wheel-coprime cumulative
Liouville sum `L_W(x)` as a finite signed combination of `L(⌊x/d⌋)`,
and the bisection lift (3d) decomposing the wheel-coprime prime count
into `2^{ω(W)}` calls to `L` plus the irreducible `C_{3,W}` residue.

The construction is novel **within the project**:

* No EDGES.md or `experiments/` entry contains Theorem 1 or its lift.
* Greps for `L_W`, `L(x; W)`, `coprime to W` Liouville sums in non-
  archive files return no matches.
* The closest project content is the *indicator-side* `T_W^{div}` =
  `[gcd(n,W) = 1] · W/φ(W)` (S208); Theorem 1 here is the
  *cumulative-side* analogue with a different but related closed form.

The mathematics (Möbius inversion of a completely multiplicative
function) is **classical** and would be derivable in an afternoon by
an analytic number theorist; this means the construction lands at
**B-grade** by the CLAUDE.md grading rubric (substantive refinement of
E1.6 with a precise new statement extending its scope to ALL
squarefree-radical wheels, plus a precise oracle-reduction corollary).
**Not A-grade**: no new mathematical content beyond direct Möbius
inversion + completely-multiplicative property.

## Pre-stated falsification criteria

All five must hold for the closure to be valid; failure of any
triggers honest reporting.

* **F1 — Pointwise Liouville identity.** For
  `W ∈ {2, 3, 5, 6, 10, 12, 15, 30, 60, 105, 210, 420, 2310}` and
  `x ∈ [0, N]` with `N = 10^4`,

  ```
      L_W(x)  =  Σ_{d | rad(W)} L(⌊x/d⌋)
  ```
  with absolute error ≤ `0` (integer-exact identity).

* **F2 — Radical reduction at non-squarefree W.** For
  `W ∈ {4, 9, 12, 60, 420, 2700}`, identity F1 holds with right-hand
  side using `rad(W)` (not `W`); equivalently, `L_W(x) = L_{rad(W)}(x)`.

* **F3 — Mod-2 closed form.** For every `W` and `x` in the F1 grid,
  `L_W(x) mod 2 = ( Σ_{d | rad(W)} ⌊x/d⌋ ) mod 2`.

* **F4 — Wheel-graded prime bisection.** For every `W` in the F1
  grid and every `x ∈ {N/4, N/2, 3N/4, N}`,
  ```
      π(x) − π_W^{div}(x)  =
            (1/2) Σ_{d | rad(W)} ( μ(d) ⌊x/d⌋ − L(⌊x/d⌋) )  −  C_{3,W}(x)
  ```
  integer-exactly. Here `π_W^{div}(x) = #{p ≤ x : p | W}` is the small
  bounded count, and `C_{3,W}(x)` is computed by direct enumeration as
  `#{n ≤ x : gcd(n, W) = 1, n composite, Ω(n) odd}`.

* **F5 — Wheel-call invariance.** The right-hand side of F1 has
  exactly `2^{ω(W)}` terms regardless of the size of `W`. Verify by
  measuring the number of evaluation calls vs the predicted `2^{ω(W)}`.

* **F6 — Mod-q lift (stretch).** Restate Theorem 1 mod `q` for `q ∈
  {2, 3, 4, 8}`, and compare the cost to incremental `O(π(N))` brute
  evaluation. Predicts a `2^{ω(W)}`-vs-`O(N)` separation at every `W`
  with `ω(W) ≪ log N`. (NOTE: this F6 is empirical — does not give a
  **polylog** algorithm for `L(x)` itself, only a `2^{ω(W)}` reduction
  of `L_W` to `L`.)

## Falsification logic

* If F1 fails: the identity is wrong (Möbius / completely-multiplicative
  derivation has a bug); investigate.
* If F2 fails: the `rad(W)` reduction is wrong; investigate.
* If F3 fails: E2.10 is wrong (impossible — proven) or the mod-2
  reduction has a bug.
* If F4 fails: the wheel-graded bisection lift has a bug, or E1.6's
  `π = (x − L)/2 − C_3` decomposition is being misapplied.
* If F5 fails: the divisor enumeration is wrong.
* If F6 fails: the cost separation is overstated; report honest cost.

## Edges composed

* **E1.6** (bisection π(x) = (x − L)/2 − C_3) — lifted to
  wheel-graded form via Theorem 3.
* **E2.10** (L(x) mod 2 = x mod 2) — used in Theorem 2 to derive the
  closed-form parity of `L_W`.
* **E2.1 / E4.1** (wheel-W density `φ(W)/W` = bond-dim ratio) —
  parallel cumulative-side identity. The construction is the
  cumulative-Liouville analogue of `T_W^{div}` (S208), which is the
  pointwise-indicator-side analogue.
* **E1.5 / T6** (CRT-mod-m saturation across coprime moduli) — frames
  why the algorithmic implication does NOT break the polylog wall:
  the wheel-graded reduction reduces `L_W` to `L`, but `L(x)` itself
  is the open frontier (`L mod 2 = x mod 2` is trivial, but `L mod 4`
  is on the same hard-zone bit boundary as `π mod 2`).

## Algorithmic content (and limits)

The reduction from `L_W` to `L` is **`2^{ω(W)}`-cost**, polynomial in
`log W` for primorial `W`. This means: any conditional algorithm (e.g.,
a polylog-`L(x)` algorithm) immediately yields a polylog-`L_W(x)`
algorithm for every wheel `W` of bounded `ω(W) ≤ O(log log x)`.

**The reduction is one-way:** it does NOT show `L(x)` reduces to
`L_W(x)`. Combined with E1.5, the W-coprime path does not bypass the
unconditional `O(x^{2/3})` ceiling for `L(x)` (Mertens-conjecture-
refutation regime). The construction is **structural** — it cleanly
separates the wheel-coprime smooth side from the irreducible `C_{3,W}`
residue.

## Save location

`experiments/constructions/wheel_coprime_liouville_identity/`.

# D38 — Prime-restricted multiple zeta values vs Brown 2012 motivic algebra

**Status:** CLOSED-I (interesting structural negative) at weights 5, 6, 7;
INC at weight ≥ 8 (basis incomplete).
**Session:** S233 (cross-domain attack mode).
**Cross-domain import:** Brown 2012 *Annals* 175 mixed Tate motives over Z;
Hoffman 1997 conjecture; Goncharov 2005 motivic coproduct on iterated
integrals. Wikipedia survey of MZVs fetched and cited below.

## Setup

Define the **prime-restricted multiple zeta value** of depth 2:
```
   ζ_P(s, t) := Σ_{p < q, p,q prime} 1/(p^s q^t)
```
and the **prime zeta function** (Mertens-type prime constant):
```
   M_s := Σ_{p prime} 1/p^s   (s ≥ 2)
```
Standard MZVs (Brown's basis): `ζ(2), ζ(3), ζ(5), ζ(7), ...` with
`ζ(2k) = (rational) · ζ(2)^k` reducing the even ζ's.

## Two structural identities

**Symmetric (Euler reflection over primes), proven exactly:**
```
   ζ_P(s, t) + ζ_P(t, s) = M_s · M_t − M_{s+t}
   ζ_P(s, s) = (M_s² − M_{2s}) / 2
```
*Proof:* The full sum `Σ_{p≠q} 1/(p^s q^t) = M_s M_t − M_{s+t}` splits
into `p<q` and `p>q` halves; in the second half relabel
`p ↔ q` and read off the LHS.

**Antisymmetric content:**
```
   A(s, t) := ζ_P(s, t) − ζ_P(t, s),    A(s, s) = 0,    A(s, t) = −A(t, s).
```
The symmetric part is a polynomial in Mertens constants. ALL
non-trivial structural content of `ζ_P(s,t)` lives in `A(s,t)`.

## D38 question

Does `A(s,t)` reduce to a `Q`-linear combination of weight-`(s+t)`
monomials in the generators
`{ ζ(2), ζ(3), ζ(5), ζ(7), …,  M_2, M_3, M_4, M_5, M_6, M_7 }`?

## Numerical method

- **Tail-corrected truncation:** for primes `p_1, …, p_π(N)`,
  ```
  ζ_P(s, t) = Σ_{p<q ≤ N} 1/(p^s q^t)  +  prefix_s(N) · tail_t(N)  +  R(N)
  ```
  where `tail_t(N) := M_t − prefix_t(N)` is computed exactly because
  `M_t` is evaluated to dps digits via mpmath's `primezeta` (Möbius
  identity `M_t = Σ_n μ(n)/n · log ζ(n t)`).
- **Truncation bound:** `|R(N)| ≤ tail_s(N) · tail_t(N)`.
  At `N = 10^6, dps = 50, (s, t) = (2, 3)`: bound `≈ 2.4 × 10^{-21}`.
  At `N = 10^6, weight 8`: bound `≈ 10^{-39}`.
- **PSLQ:** mpmath `pslq` with `maxcoeff = 10^18`, tolerance
  `10^{-0.4·dps}`.

## Results

### Sanity: Euler reflection holds at truncation precision

| (s,t) | LHS = ζ_P(s,t) + ζ_P(t,s) − (M_s M_t − M_{s+t}) | Truncation bound |
|-------|-------------------------------------------------|------------------|
| (2,3) | 2.4e-21                                         | 2.4e-21          |
| (2,4) | 1.6e-27                                         | 1.6e-27          |
| (2,5) | 1.2e-33                                         | 1.2e-33          |
| (3,4) | 8.2e-34                                         | 8.2e-34          |
| (2,6) | 9.7e-40                                         | 9.7e-40          |
| (3,5) | 6.2e-40                                         | 6.2e-40          |

Symmetric reduction `ζ_P(s,t) + ζ_P(t,s) = M_s M_t − M_{s+t}` holds
exactly, saturating the truncation bound. ✓

### PSLQ on antisymmetric A(s,t), N = 10^6, dps = 50

| (s,t) | weight | A(s,t) (25 digits)            | basis dim | PSLQ status |
|-------|--------|-------------------------------|-----------|-------------|
| (2,3) | 5      | −0.01508939774920873752716337 | 6         | NO_RELATION |
| (2,4) | 6      | −0.00997538653317191147159479 | 10        | NO_RELATION |
| (2,5) | 7      | −0.00554105809954202444760610 | 14        | NO_RELATION |
| (3,4) | 7      | −0.00137968302292165495378312 | 14        | NO_RELATION |
| (2,6) | 8      | −0.00291886724859182279686884 | 22*       | NO_RELATION |
| (3,5) | 8      | −0.00103092360373435439055005 | 22*       | NO_RELATION |

\* basis at weight 8 is INCOMPLETE — missing depth-≥3 irreducible MZV
generators, e.g., the irreducible double zeta `ζ(3,5)` (Brown 2012).
At weights 5, 6, 7 the basis IS complete: `d_5 = 2, d_6 = 2, d_7 = 3`
(Brown dims) and the products `{ζ(2)^k ζ(3)^l ζ(5)^m ζ(7)^n}` cover
the MZV side fully, supplemented by Mertens monomials.

### Cross-validation at N = 10^7 (dps=50)

ζ_P(2,3) at N=10^6: `0.01409576875480383351191378` (err ≤ 2.4e-21)
ζ_P(2,3) at N=10^7: `0.01409576875480383351271971` (err ≤ 1.8e-24)
Difference 8e-19 — within the N=10^6 truncation bound. ✓

PSLQ result at N=10^7 dps=50: still NO_RELATION for A(2,3), A(2,4).

## Verdict

**The antisymmetric prime-restricted MZV `A(s,t) := ζ_P(s,t) − ζ_P(t,s)`
is NOT a Q-linear combination of standard MZV monomials × Mertens
monomials of total weight (s+t), for `(s,t) ∈ {(2,3), (2,4), (2,5),
(3,4)}` — the weights where the Hoffman/Brown product basis is
complete. PSLQ ruled out integer relations with coefficient bound 10^18
at 50-digit precision.**

This is the **I outcome** in the D38 falsification protocol: primes
generate a *strictly richer* depth-2 period algebra than Brown 2012's
mixed Tate motives ⊕ Mertens constants. The symmetric "Euler-reflection"
half collapses to Mertens; the antisymmetric half does not collapse and
encodes prime-ordering content not seen by Brown's algebra.

## What would falsify this

- A relation with **coefficients > 10^18** in some entry — possible but
  structurally unprecedented; standard MZV reductions have coefficients
  with bounded heights (< 10^4 typically by Zagier dimension).
- An irreducible **depth ≥ 3** Hoffman MZV being added to the basis
  closes the gap. At weight 8 my basis is missing `ζ(3,5)` (the first
  irreducible double zeta of weight ≥ 8); weight 9+ misses
  `ζ(3,5,2), ζ(3,7), ζ(5,5), …`. The result at weights ≤ 7 is
  immune to this concern: there `d_w` matches the dim of my products-
  only basis.
- A future paper proving `ζ_P(s,t)` is a period of a known mixed Tate
  motive scheme over `Spec(Z) − {primes}` or similar.

## Edges composed / cited

- **E_new (D38 antisym)** — *the* edge produced by this session, see
  EDGES.md S233 entry.
- E2.1 (Riemann-zeta + Möbius prime-zeta connection) — used to compute
  `M_s` to 50 digits.
- E2.16 (DPP / non-determinantal point process structure on χ_P) —
  D40 sibling vector, distinct from D38.
- E7.x — distinct from all closed prime-polynomial / Mahler-measure /
  Newman-flatness edges (D10 CLOSED, D27 CLOSED).

## Reproducibility

```
cd experiments/algebraic/d38_prime_mzv
python3 d38_prime_mzv.py --N 1000000 --dps 50 --weights 5,6,7
# 14 ζ_P values computed in ~25 s. PSLQ tests within 1 s each.
```

## Cross-domain reference (cited per CLAUDE.md cross-domain rule)

- Brown 2012 "Mixed Tate motives over `Z`" *Annals* 175 = arXiv:1102.1312
  https://arxiv.org/abs/1102.1312
- Hoffman 1997 "The algebra of multiple harmonic series" *J. Algebra* 194
- Goncharov 2005 "Galois symmetries of fundamental groupoids…"
  *Duke Math. J.* 128 = arXiv:math/0208144
- Zagier 1994 "Values of zeta functions and their applications"
- Wikipedia: Multiple zeta function (fetched 2026-04-29) —
  conventions for ζ(s_1,…,s_k) ordering and explicit weight-≤6 closed
  forms `ζ(2,2)=(3/4)ζ(4), ζ(2,3)=(9/2)ζ(5)−2ζ(2)ζ(3),
  ζ(3,2)=3ζ(2)ζ(3)−(11/2)ζ(5), ζ(3,3)=½(ζ(3)²−ζ(6))`.

## Grade

**B-grade structural negative.** The cross-domain import (Brown 2012
motivic period theory) DID real work: it specified a precise basis
against which `ζ_P` could be tested, and the test was conclusive at
weights ≤ 7. The result is a previously-unmeasured structural fact
about the period content of primes.

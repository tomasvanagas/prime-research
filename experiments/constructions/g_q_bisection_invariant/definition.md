# Construction: g_q — q-generalised Liouville bisection

**Composition target:** C1 from `NOVELTY_CHALLENGES.md` §1
**Edges composed:** **E1.6** (Liouville bisection of pi(x) mod 2 at q=2)
                    **E1.5** (per-step conditional entropy saturates at
                              h_2(pi(X)/X), sharpened in S69 for pi(x) mod m)

## Signature

For q >= 2 a positive integer, define the *paired* counting map

```
   g_q : [1, N] -> Z/qZ x Z/qZ
   g_q(x) = ( A(x) mod q, C_3(x) mod q )
```

where

```
   A(x)   = #{ n <= x : Omega(n) odd }            = (x - L(x)) / 2
   C_3(x) = #{ n <= x : Omega(n) odd, >= 3 }      = A(x) - pi(x)
   L(x)   = sum_{n <= x} lambda(n)                (Liouville summatory)
```

so that pi(x) = A(x) - C_3(x) holds in Z (E1.6, verified bit-exact for
x in [1, 2 * 10^6]).

## Why this is a composition (not a duplicate)

E1.6 records the q=2 case of g_q as a pair (A mod 2, C_3 mod 2) with three
properties:
  (i)  each component uniform on {0, 1};
  (ii) mutual info I(A mod 2 ; C_3 mod 2) = 1.94 * 10^-5 bits (independent);
  (iii) XOR is pi(x) mod 2.

E1.5 (S69-sharpened) records the per-step conditional entropy saturation
H( pi(x) mod m | pi(x-1) mod m ) = h_2( pi(X) / X ) + O( 1 / pi(X) ),
in regime m << pi(X), uniformly in m.

The composition's content. The same closed-form mechanism that drives E1.5
(monotone counter, increment in {0, 1}, conditional independence of the
indicator from the state) applies to *any* monotone integer-valued counter.
Applied to A(x) and C_3(x) it predicts:

```
   H( A(x)   mod q | A(x-1)   mod q )  =  h_2( A(X) / X )    + O(1/A(X))
   H( C_3(x) mod q | C_3(x-1) mod q )  =  h_2( C_3(X) / X )  + O(1/C_3(X))
   H( pi(x)  mod q | pi(x-1)  mod q )  =  h_2( pi(X) / X )   + O(1/pi(X))
```

But unlike pi, the densities A(X)/X and C_3(X)/X **do not vanish**:

```
   A(X)   / X  ->  1/2          (Pillai-Selberg, density of Omega odd)
   C_3(X) / X  =   A(X)/X - pi(X)/X  ->  1/2     (since pi(X)/X -> 0)
```

So the prediction is

```
   lim_{X -> inf}  H( A   mod q | prev )  =  h_2(1/2)  =  1 bit / step
   lim_{X -> inf}  H( C_3 mod q | prev )  =  h_2(1/2)  =  1 bit / step
   lim_{X -> inf}  H( pi  mod q | prev )  =                0 bits / step
```

Each component of the bisection carries asymptotically **one full bit per step**,
yet their integer difference pi = A - C_3 carries asymptotically **zero
bits per step**. Whether E1.5's q-independence (its constancy across moduli
in regime m << density) generalises to A and C_3 is the empirical question.

## Joint per-step rate prediction

Per-step increments (delta_A, delta_C_3) take values in three (not four)
states because delta_C_3 = 1 forces delta_A = 1 (incompatibility of
"Omega odd >= 3" with "Omega = 1"):

```
   (delta_A, delta_C_3) = (0, 0) iff Omega(x) even   prob ~ 1/2
                         = (1, 0) iff Omega(x) = 1   prob = pi(X)/X
                         = (1, 1) iff Omega(x) odd >= 3   prob ~ 1/2 - pi(X)/X
```

so the joint per-step conditional entropy of g_q is, in the regime q^2 << pi(X)
(state space coverage),

```
   H( g_q(x) | g_q(x-1) )  =  H_3( 1 - A(X)/X,  pi(X) / X,  C_3(X) / X )
                             + O( 1 / pi(X) )
```

This is the construction's *novel quantitative prediction*. It is a
strict generalisation of E1.5 (recovered as the q=1 marginal of the
A-C_3 component) and a strict generalisation of E1.6 (q=2 case as the
two-state marginal).

## Marginal q-stable independence prediction (E1.6 generalisation)

E1.6 records I(A mod 2 ; C_3 mod 2) ~ 0 over x ~ Uniform[1, X]. The
construction predicts the analogue for q in {3, 5, 7, 11, 13}:

```
   I( A mod q ; C_3 mod q )  -> small (< 0.01 bits at X = 2 * 10^6)
```

If true for all q tested, that is a q-stable strengthening of E1.6.
If false for some q > 2, the failure is genuinely new structure
(non-binary correlation between the two summatory parities).

## What the construction is NOT

* It is not a new pseudorandomness measure — these are exact closed-form
  predictions, not floor-of-noise tests. The output is a verdict on the
  predictions, not a residual.
* It is not a path to faster pi(x) — knowing g_q gives pi mod q via
  subtraction, but g_q is no easier than pi to compute (each component
  requires Omega-stratified counting).
* It is not a duplicate of E1.6 — E1.6 is q=2 only.
* It is not a duplicate of E1.5 — E1.5 is for pi(x) mod m, not for the
  bisection components A, C_3.

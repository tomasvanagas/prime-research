# T_Q triple correlation = truncated Hardy-Littlewood prime triple
# singular series (S209 construction, three-point companion to S205)

## Object defined

For positive integer N, integer Q >= 1 and shifts h_1, h_2 (integers), the
**triple T_Q correlation** is

```
    R_{h_1, h_2}(Q, N)        :=  <T_Q(n) T_Q(n+h_1) T_Q(n+h_2)>_n,
    R_{h_1, h_2}^{conn}(Q,N)  :=  R_{h_1,h_2}
                                   - 3 disconnected-pair pieces
                                   + 2 mu^3
                                   (subtract Wick contractions to isolate
                                   the genuine triple cumulant)
                                = K3(Q, h_1, h_2),
```

with `mu = <T_Q> = pi(N)/N + O(1/N)` (the rank-1 wheel mode).

## Closed-form prediction

Two parallel statements, each verified empirically and derived
analytically below.

### A. Divisor-restricted (primorial conductor W) - EXACT

Let `T_W^{div}(n) := (pi(N)/N) * Sum_{q | W, q sqf} mu(q)/phi(q) * c_q(n)`.
By S208's Mobius collapse this equals `(pi(N)/N) * (W/phi(W)) *
[gcd(n, W) = 1]` for any squarefree W. For W primorial, `T_W^{div} = T_W`
(the divisors of primorial W ARE the squarefree integers <= W with all
prime factors among the primes of W).

**Theorem (S209, exact at primorial W; pointwise inclusion-exclusion).**
For every primorial W >= 1 and integers (h_1, h_2),

```
    <T_W^{div}(n) T_W^{div}(n+h_1) T_W^{div}(n+h_2)>_n
        =  (pi(N)/N)^3  *  prod_{p | W} (p - nu_p(0, h_1, h_2)) p^2 / (p-1)^3
        =  (pi(N)/N)^3  *  S_HL^{(W)}(0, h_1, h_2)
```

where `nu_p(0, h_1, h_2) := #{distinct residues among {0, h_1, h_2} mod p}`
and `S_HL^{(W)}` is the **conductor-W truncation of the Hardy-Littlewood
prime triple singular series**. The full (untruncated) HL prediction is
`S_HL(0, h_1, h_2) = prod_{p prime} (1 - nu_p / p) (1 - 1/p)^{-3}`.

### B. General Q (squarefree-cap) - APPROXIMATE

Let `T_Q(n) := (pi(N)/N) * Sum_{q sqf <= Q} mu(q)/phi(q) * c_q(n)`. By
Ramanujan-Fourier orthogonality (prime-by-prime factorization for the
diagonal piece plus 2-point cross-pairs and a trivariate cross-pair),
the triple correlation factorises asymptotically as

```
    R_{h_1,h_2}(Q, N)
        =  (pi(N)/N)^3  *  prod_{p <= Q (squarefree)} G_p(h_1, h_2)
           +  O(N^{-1+eps})
```

where the prime-p contribution is

```
    G_p(h_1, h_2)  =  1
                    +  (1/(p-1)^2) * [ c_p(h_1) + c_p(h_2) + c_p(h_2 - h_1) ]
                    +  (-1/(p-1)^3) * f_p(h_1, h_2)
```

with `f_p(h_1, h_2) = (1/p) Sum_{r mod p} c_p(r) c_p(r+h_1) c_p(r+h_2)`
given by the **closed form**:

```
    f_p =  (p-1)(p-2)   if nu_p({0, h_1, h_2}) = 1   (all coincide mod p)
        =  -(p-2)       if nu_p = 2                  (one pair coincides)
        =  +2           if nu_p = 3                  (all distinct mod p)
```

**Identity.** A direct algebraic check (proven below in `proof.md`-style
in tq_triple_correlation_results.md) shows

```
    G_p(h_1, h_2)  =  (p - nu_p(h_1, h_2)) * p^2 / (p-1)^3
                  =  (1 - nu_p/p) * (1 - 1/p)^{-3}
```

i.e. **the prime-p factor of G_p IS the prime-p factor of the
Hardy-Littlewood prime triple singular series, exactly**.

## Reduction to lower-order identities

* **h_2 = 0** (degenerate triple to pair `{0, h_1}`): `nu_p(0, h_1, 0) =
  nu_p(0, h_1) = 1` if `p | h_1` else `2`. The triple G_p collapses
  algebraically to (W/phi(W)) * (the pair G_p) since
  T_W^{div}(n)^2 = (W/phi(W)) * (pi/N) * T_W^{div}(n). Cross-checked
  empirically.
* **h_1 = h_2 = 0**: `<T^3>` = (W/phi(W))^2 (pi/N)^2 * <T> =
  (W/phi(W))^2 * (pi/N)^3 (verified F4-style at h_1 = h_2 = 0).

## Pre-stated falsification criteria

| Tag | Description | Pass band |
|-----|-------------|-----------|
| F1 | Identity ratio R_{h_1,h_2}(Q, N) / [(pi/N)^3 prod_{p|W} G_p(h_1, h_2)] for primorial W | within 0.5% on every (W, h_1, h_2) at d >= 18 |
| F2 | Identity ratio R_{h_1,h_2}(Q, N) for general Q vs prod_{p sqf <= Q} G_p | within 1.5% at d=20, Q in {30, 210, 2310, sqrt(N)} |
| F3 | HL recovery: ratio of R_{h_1,h_2}(Q=sqrt(N), N) to (pi/N)^3 S_HL(0,h_1,h_2) for admissible (h_1, h_2) | within 1% at d=20 |
| F4 | h_1 = h_2 = 0 self-consistency: <T_W^3> = (pi/N)^3 (W/phi(W))^2 | within 0.01% (this is a triviality check) |
| F5 | Reduction to S205 2-point at h_2 = 0: prediction comparison with S205 | within 0.01% (algebraic identity) |
| F6 | Inadmissible triple (some prime p has nu_p = p, e.g. (h_1, h_2) = (2, 4) and p=2: nu_2({0,2,4}) = 1) | exact zero from prime-density side; T_Q triple is **non-zero** but matches G_p formula |

## Failure modes (negative cases)

If the construction's three-point function does not match the
HL prime triple singular series at any (W, h_1, h_2), the identity
collapses; this would imply the S205 two-point identity is special-case
and not part of a coherent k-point family.

If G_p does not factor cleanly via the case-by-case f_p table, the
Ramanujan-Fourier derivation has a bug.

## What this construction adds beyond S205 + S208

* **Three-point structure.** S205 covered the 2-point case (twin-prime
  HL singular series). S209 extends to the 3-point case (prime triple
  HL singular series), giving a closed form at primorial conductors and
  prime-by-prime factorization at general Q.

* **Composes E2.16's "DPP failure" 3-point fact with E2.13's HL match.**
  E2.16 noted the 3-point HL singular series factors over PRIMES (not
  pairs) — ruling out kernel-DPP/PPP descriptions. S209 exhibits an
  explicit pointwise scalar field whose 3-point function equals exactly
  this prime-factored HL series.

* **Frees the wheel-W -> HL connection from purely-pair statements.**
  E2.1 / E4.1 wheel-W gain is per-pair; S209 shows the wheel indicator's
  3-point function realizes the prime-3-tuple HL series at conductor
  truncation.

## Files

* `tq_triple_correlation.py` — main script.
* `tq_triple_correlation_results.md` — empirical verification + analytic proof.
* `tq_triple_correlation_results.json` — raw numerics.

## Edges this construction touches

* **Refines E2.1** (T_Q spike content) with the **3-point** form of the
  spike subspace's correlation.
* **Extends E2.13** (Gowers U^k -> HL singular series) to the
  **prime-triple** singular series, with explicit pointwise object T_W^{div}.
* **Composes with E2.16** (DPP failure: 3-point HL factors over primes)
  by providing the explicit factor-by-factor identity.
* **Refines C9 / S191 / S205 / S208** by extending the spike-pointwise
  framework to one higher correlation order.
* **Uses E1.6 / E2.2** (parity bisection) at the q=2 prime factor.

## NO cross-domain technique imported

This construction uses only project-internal mathematics: Ramanujan
sums, Mobius inversion, Hardy-Littlewood singular series (already
catalogued in EDGES.md), and direct enumeration. No new named technique
is added to CROSS_DOMAIN_TECHNIQUES.md.

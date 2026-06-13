# T_Q quadruple correlation = truncated Hardy-Littlewood prime quadruple
# singular series (S?, four-point companion to S205 / S208 / S209)

## Object defined

For positive integer N, integer Q >= 1 and shifts h_1, h_2, h_3, the
**quadruple T_Q correlation** is

```
    R_{h_1, h_2, h_3}(Q, N)
        :=  <T_Q(n) T_Q(n+h_1) T_Q(n+h_2) T_Q(n+h_3)>_n.
```

with the Ramanujan-spike pointwise approximator

```
    T_Q(n) := (pi(N)/N) Sum_{q sqf <= Q} (mu(q) / phi(q)) c_q(n)
```

(S191 / C9). At primorial conductors W, the divisor-restricted variant
is

```
    T_W^{div}(n) := (pi(N)/N) * (W / phi(W)) * [gcd(n, W) = 1]
```

(S208 collapse).

## Closed-form prediction

### A. Divisor-restricted (primorial conductor W) - EXACT pointwise identity

**Theorem (this construction, exact at squarefree W).**
Let W >= 1 be squarefree.  For all integers (h_1, h_2, h_3) and all N
sufficiently large,

```
    <T_W^{div}(n) T_W^{div}(n+h_1) T_W^{div}(n+h_2) T_W^{div}(n+h_3)>_n
        =  (pi(N)/N)^4
           * prod_{p | W} (p - nu_p(0, h_1, h_2, h_3)) * p^3 / (p-1)^4
        =  (pi(N)/N)^4 * S_HL^{(W)}(0, h_1, h_2, h_3),
```

where `nu_p(0, h_1, h_2, h_3) := #{distinct residues among
{0, h_1, h_2, h_3} mod p}` and `S_HL^{(W)}` is the **conductor-W
truncation of the Hardy-Littlewood prime quadruple singular series**.

For non-squarefree W, the same identity holds with W replaced by its
radical `rad(W)`.

**Proof.** Direct from S208's pointwise identity:

  T_W^{div}(n) = (pi(N)/N) * (W/phi(W)) * [gcd(n, W) = 1].

Hence

  prod_{i=0..3} T_W^{div}(n+h_i)
    = (pi(N)/N)^4 * (W/phi(W))^4
      * prod_{i=0..3} [gcd(n+h_i, W) = 1].

The product of indicators is multiplicative across the primes p | W:

  prod_{i=0..3} [gcd(n+h_i, W) = 1]
    = prod_{p | W} prod_{i=0..3} [n+h_i not = 0 mod p]
    = prod_{p | W} [n not in -{0, h_1, h_2, h_3} mod p].

For each prime p | W, the density (over residues r mod p) of "r not in
-{0, h_1, h_2, h_3} mod p" is exactly `(p - nu_p) / p`.  Hence the
n-average factorises as

  <prod_{i=0..3} [gcd(n+h_i, W) = 1]>_n
    = prod_{p | W} (p - nu_p) / p
    + O(1/N).

Substituting `(W/phi(W))^4 = prod_{p | W} (p/(p-1))^4`,

  (W/phi(W))^4 * prod_{p|W} (p - nu_p)/p
    = prod_{p|W} (p - nu_p) * p^3 / (p-1)^4.

That is the prime-p factor of the Hardy-Littlewood prime *quadruple*
singular series.  Q.E.D.

### B. General Q (squarefree-cap) - APPROXIMATE prime-by-prime form

For general squarefree-capped Q, the diagonal Ramanujan-Fourier mode
(all four conductors equal a single prime p) gives the same prime-p
factor

```
    G_p^{(4)}(h_1, h_2, h_3)  =  (p - nu_p) * p^3 / (p-1)^4.
```

Off-diagonal terms (mixed conductors) contribute O(N^{-1+eps}) finite-N
errors.  The expected ratio

```
    R_{h_1,h_2,h_3}(Q, N) / [(pi/N)^4 * prod_{p sqf primes <= Q} G_p^{(4)}]
        -> 1 + O(Q*log Q / N).
```

### Ramanujan-Fourier 4-point cross-check (G_p^{(4)} algebraic identity)

For sanity we also implement the direct expansion over the
4-cumulant of T_p centred about its rank-1 mean:

```
    G_p^{(4)}(h_1, h_2, h_3)
        =  1
        +  1/(p-1)^2 * Sum_{i<j in {0,h_1,h_2,h_3}} c_p(j - i)
        -  1/(p-1)^3 * Sum_{|S|=3 \subset {0,h_1,h_2,h_3}} f_p^{(3)}(S)
        +  1/(p-1)^4 * f_p^{(4)}(h_1, h_2, h_3),
```

where `f_p^{(k)}(H)` is the direct sum `(1/p) sum_{r mod p}
prod_{h in H} c_p(r + h)`.  This expansion must equal
`(p - nu_p) p^3 / (p-1)^4` for every prime p and every triple
(h_1, h_2, h_3) (F7 sanity).

## Reduction to lower-order identities

* **h_3 = 0**: `T_W^{div}(n)^2 = (W/phi(W)) * (pi/N) * T_W^{div}(n)`,
  so the 4-point correlation collapses to `(W/phi(W)) * <T^3>_h_1,h_2`
  (S209 prediction multiplied by an overall constant).  Cells (2,4,0)
  / (2,6,0) / (6,12,0) verify F5.

* **h_2 = h_3 = 0**: collapses to `(W/phi(W))^2 * <T^2>_h_1`
  (S205 prediction multiplied by an overall constant^2).  Cells
  (2,0,0) / (6,0,0) verify F5.

* **h_1 = h_2 = h_3 = 0**: `<T_W^{div^4}> = (pi/N)^4 *
  (W/phi(W))^3` (per S208's pointwise scalar field);
  the closed-form prediction
  `(pi/N)^4 * prod_{p|W} (p-nu_p) p^3/(p-1)^4 = (pi/N)^4 *
  prod_{p|W} (p-1) p^3/(p-1)^4 = (pi/N)^4 * prod_{p|W} p^3/(p-1)^3 =
  (pi/N)^4 * (W/phi(W))^3`.  F4 algebraic identity.

## Pre-stated falsification criteria

| Tag | Description | Pass band |
|-----|-------------|-----------|
| F1  | Identity ratio R_{h1,h2,h3}(W, N) / [(pi/N)^4 prod_{p|W} G_p^{(4)}] for primorial W. | within 0.5% on every (W, h_1, h_2, h_3) at d >= 18. |
| F2  | Identity ratio R_{h1,h2,h3}(Q, N) for general Q vs prod_{p sqf <= Q} G_p^{(4)}. | within 2.5% at d=20, Q in {30, 210, 2310, sqrt(N)}. |
| F3  | HL recovery: ratio of R_{h1,h2,h3}(Q=sqrt(N), N) to (pi/N)^4 S_HL(0,h_1,h_2,h_3) for admissible (h_1,h_2,h_3). | within 3% at d=20. |
| F4  | h_1 = h_2 = h_3 = 0 self-consistency: <T_W^4> = (pi/N)^4 (W/phi(W))^3. | within 0.01% (algebraic). |
| F5  | Reduction to S209 3-point at h_3 = 0: <T_W^{div^4}_h_1,h_2,0> = (W/phi(W)) * <T_W^{div^3}_h_1,h_2>. | within 0.01% (algebraic). |
| F6  | Inadmissible quadruple at p with p|W: closed-form prediction = 0; empirical primorial-W ratio = 0/0 (undefined) or empirical correlation < 1e-6. | empirical correlation magnitude bounded. |
| F7  | G_p^{(4)} closed-form algebraic consistency vs the direct Ramanujan-Fourier 4-cumulant expansion on (p, h_1, h_2, h_3) cells with p in {2,3,5,7,11,13}. | abs_err < 1e-10 on every cell. |

## Failure modes

If F1 fails for primorial W: the primorial-W S208 pointwise collapse is
incompatible with the 4-point HL singular series — either the proof is
wrong or finite-N corrections at d=18..20 are larger than predicted.
Either outcome is genuine new content.

If F7 fails: the Ramanujan-Fourier 4-cumulant expansion is being
mis-derived; the closed form `(p - nu_p) p^3 / (p-1)^4` is independently
provable via S208's pointwise route, so an F7 failure is a bug in the
expansion, not a flaw in the closed form.

If F3 fails dramatically: the spike-pointwise framework's correlation
hierarchy diverges from the HL prime-tuple singular series at the
4-tuple level, suggesting the framework breaks down for k >= 4 — an
A-grade-shape negative result.

## What this construction adds beyond S205 + S208 + S209

* **Four-point structure.** S205 covered the 2-point case (twin-prime
  HL); S208 the 1-point pointwise indicator; S209 the 3-point case
  (prime triple HL).  This construction extends to the 4-point case
  (prime quadruple HL singular series), giving an explicit closed form
  at primorial conductors and a direct generalisation of the
  prime-by-prime factorisation to one higher correlation order.

* **Composes E2.16's 3-point HL prime factorisation with the cube-shift
  case of E2.13 (Gowers U^2).**  The 4-point HL singular series for
  shifts (0, h_1, h_2, h_1 + h_2) is the U^2 norm of T_Q on the
  cube vertex set; this construction subsumes that case.

* **Establishes the k-point spike-pointwise -> HL identity hierarchy
  through k=4.**  Combined with the now-proven k = 1, 2, 3 statements
  (S208 / S205 / S209), this gives the first verified four-step
  hierarchy of pointwise spike correlation identities matching the
  Hardy-Littlewood prime k-tuple singular series at primorial
  conductors.

## Files

* `tq_quad_correlation.py` — main script.
* `tq_quad_correlation_results.md` — empirical verification + analytic notes.
* `tq_quad_correlation_results.json` — raw numerics.

## Edges this construction touches

* **Refines E2.1** (T_Q spike content) with the **4-point** form of the
  spike subspace's correlation.
* **Extends E2.13** (Gowers U^k -> HL singular series) to the
  **prime-quadruple** singular series, with explicit pointwise object
  T_W^{div}.
* **Composes with E2.16** (DPP failure: 3-point HL factors over primes)
  by extending the prime-by-prime factor identity to 4-tuples.
* **Composes with E1.6 / E2.2** at the q=2 prime factor (parity
  bisection): all admissible 4-tuples have all h_i even.
* **Refines C9 / S191 / S205 / S208 / S209** by extending the
  spike-pointwise framework to one higher correlation order.

## NO cross-domain technique imported

This construction uses only project-internal mathematics: Ramanujan
sums, Mobius inversion, Hardy-Littlewood singular series (already
catalogued in EDGES.md), and direct enumeration.  No new named
technique is added to CROSS_DOMAIN_TECHNIQUES.md.

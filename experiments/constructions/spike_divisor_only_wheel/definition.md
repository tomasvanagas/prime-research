# Divisor-only restriction of the Mertens-Ramanujan spike intensity

## Object

Let `W ≥ 2` be a positive integer. Define the **divisor-only spike intensity**
at conductor cutoff `W`,

```
    M_W^{div}(n)  :=  Σ_{q | W, q squarefree}  mu(gcd(q, n)) / phi(q / gcd(q, n))
```

and the corresponding pointwise **divisor-only spike approximator** of `chi_P`,

```
    T_W^{div}(n)  :=  (pi(N) / N) · M_W^{div}(n).
```

This is the C9 / S191 construction (`ramanujan_spike_pointwise/`) restricted
to summing only over **divisors of `W` that are squarefree**, instead of over
all squarefree integers `≤ W`. At general `W` the two sets differ; at primorial
`W = p_1 · p_2 · ... · p_k` they differ by exactly the squarefree integers
`≤ W` that fail to divide `W` (e.g. for `W = 6`, the index `q = 5`).

## Pre-stated identity (proven analytically — to be verified pointwise)

**Theorem (closed form for `M_W^{div}` at squarefree `W`).** For every
squarefree integer `W ≥ 1` and every positive integer `n`,

```
    M_W^{div}(n)  =  [gcd(n, W) = 1] · W / phi(W).
```

**Proof.** Write `W = p_1 · p_2 · ... · p_k`. Every squarefree divisor `q | W`
is `q_S = ∏_{p ∈ S} p` for a unique `S ⊆ {p_1, ..., p_k}`. Let
`P_n = {p ∈ {p_1, ..., p_k} : p | n}` and `P_{!n} = {p_1,...,p_k} \ P_n`.
Decompose `S = A ∪ B` with `A = S ∩ P_n` (the primes in `S` dividing `n`)
and `B = S ∩ P_{!n}`. Then `gcd(q_S, n) = ∏_A` and `q_S / gcd(q_S, n) = ∏_B`.
Hence

```
    mu(gcd(q_S, n)) / phi(q_S / gcd(q_S, n))
        = (-1)^|A| / ∏_{p ∈ B} (p - 1).
```

Summing over `S = A ∪ B` factorises:

```
    M_W^{div}(n) = ( Σ_{A ⊆ P_n} (-1)^|A| ) · ( Σ_{B ⊆ P_{!n}} ∏_{p ∈ B} 1/(p-1) ).
```

The second factor equals `∏_{p ∈ P_{!n}} (1 + 1/(p-1)) = ∏_{p ∈ P_{!n}} p/(p-1)`.
The first factor equals `(1 - 1)^|P_n| = 0` if `P_n ≠ ∅` and `1` if `P_n = ∅`.
Therefore

```
    M_W^{div}(n) = [P_n = ∅] · ∏_{p | W} p/(p - 1)
                 = [gcd(n, W) = 1] · W / phi(W).                              □
```

For `W` not squarefree, the sum only sees the **radical** `rad(W)`:

```
    M_W^{div}(n)  =  [gcd(n, rad(W)) = 1] · rad(W) / phi(rad(W)).
```

## Why this is novel relative to the parent construction (S191 / C9)

S191 derived `T_Q(n) | gcd(n, W) = 1 = pi(N)/N · W/phi(W)` for **primorial
`Q = W`**, restricted to coprime `n`, as a corollary of the full `M_Q`
sum (over all squarefree `q ≤ Q`). The divisor-only restriction
`M_W^{div}` admits a **strictly stronger statement**:

* Holds for **all squarefree `W`**, not just primorials.
* Gives the value at `non-coprime n` exactly: identically zero.
* Is **closed-form on the full domain** (the original `T_Q | gcd=1`
  identity is a single-residue-class statement).
* Reduces in a controlled way for non-squarefree `W` (depends only on
  `rad(W)`).

The numerical content of the identity is that the divisor-only restriction
collapses `T_W^{div}` to a **scaled wheel indicator** — *exactly* the
prime-density-weighted wheel sieve admissibility function. The non-divisor
squarefree terms in the full `M_W` are precisely what makes `T_W` non-flat
on the coprime cosets and capable of separating primes from composites
within those cosets.

## Pre-stated falsification criteria

Run on `N = 10^5` and `N = 10^6` (full chi_P available by sieve). All
six criteria must hold for the closure to be valid; failure of any
criterion triggers honest reporting of the failure mode.

* **F1 — Squarefree pointwise identity.** For
  `W ∈ {2, 6, 30, 210, 2310}` (primorials) and every `n ∈ [1, N]`:
  ```
      M_W^{div}(n)  =  [gcd(n, W) = 1] · W / phi(W)
  ```
  with absolute error ≤ `10^{-12}`.

* **F2 — Squarefree-non-primorial pointwise identity.** For
  `W ∈ {15, 105}` (squarefree, non-primorial: 3·5, 3·5·7):
  same identity, abs err ≤ `10^{-12}`.

* **F3 — Radical reduction at non-squarefree `W`.** For
  `W ∈ {12, 60, 420}` (with `rad(W) ∈ {6, 30, 210}`):
  ```
      M_W^{div}(n)  =  [gcd(n, rad(W)) = 1] · rad(W) / phi(rad(W))
  ```
  abs err ≤ `10^{-12}`.

* **F4 — Mean unit normalisation.** `⟨M_W^{div}⟩_n  =  1` for all
  squarefree `W` (as `(phi(W)/W) · (W/phi(W)) = 1`); for non-squarefree
  `W`, `⟨M_W^{div}⟩_n = 1` (replacing `W` by `rad(W)`).

* **F5 — L² mean-square deviation closed form.** For squarefree `W`,
  ```
      || T_W^{div} - chi_P ||² / N  =  rho · (1 - rho / phi_w)
  ```
  where `rho = pi(N)/N` and `phi_w = phi(W)/W`, up to `O(1/N)` finite-
  size error and the `O(W)` boundary term from the few primes
  `p ≤ W` that themselves divide `W`. Specifically: at `W = 210`,
  predicted ratio `rho · (1 - rho · W/phi(W))` should match measured
  `||T_W^{div} - chi_P||²/N` to within `10^{-3}` relative.

* **F6 — Comparison with full `T_W` (S191).** The full `T_W(n)`
  (`ramanujan_spike_pointwise/`) is **strictly more discriminating**
  on coprime cosets:
  ```
      Pearson( chi_P, T_W ) > Pearson( chi_P, T_W^{div} )
  ```
  for every `W ∈ {6, 30, 210, 2310}` at `N = 10^6`. (T_W^{div} on
  coprime cosets is constant, giving zero pointwise discrimination
  among coprime n.)

## Falsification logic

* If F1, F2, F3 hold: the closed-form theorem is verified pointwise.
* If F4 holds: the wheel-density normalisation is exact.
* If F5 holds: the L² behaviour matches the Bernoulli-on-coprime model.
* If F6 holds: the divisor restriction is *informationally weaker* than
  the full `T_W` — quantifying what the non-divisor squarefree
  conductors contribute.
* If any F_i fails: the closed form is wrong / the implementation has
  a bug / a deeper structure exists; investigate and report honestly.

## Edges composed

* **E2.1** (MPS bond-dim identity) — refines the wheel-density factor
  `phi(W)/W` from a rank-bound to a *pointwise* identity.
* **C9 / S191** (full `T_Q` pointwise identity) — `M_W^{div}` is the
  divisor-only restriction with strictly cleaner closed form.
* **E2.2** (Liouville/parity identity) — `M_2^{div}(n) = 2 · [n odd]`
  exactly recovers the parity bisection bit.

## Algorithmic content

* **Cost** per evaluation: `O(2^{omega(W)})` where `omega(W)` is the
  number of distinct prime divisors of `W`. For `W = 2310 = 2·3·5·7·11`
  that's `2^5 = 32` operations; far cheaper than the full `M_W` cost
  `O(W · omega(n))` for the parent construction.
* **Output**: `T_W^{div}(n) = (pi(N)/N) · W/phi(W)` if `gcd(n, W) = 1`,
  else `0`. This is just the **scaled wheel indicator** — confirming
  the divisor restriction collapses to a finite-state object with no
  internal structure on the coprime cosets.
* **No polylog opening**: the divisor-only restriction is informationally
  weaker than the full `T_Q`. The composition is mathematically exact
  but algorithmically inert (the value is determined by `gcd(n, W)`,
  computable trivially). The novel content is the **closed-form
  identity** itself, which clarifies what `M_Q` does *beyond* the
  divisor-only part.

## Save location

`experiments/constructions/spike_divisor_only_wheel/`.

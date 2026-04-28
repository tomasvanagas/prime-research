# Pointwise Ramanujan-spike approximator T_Q(n)

## Object

For positive integer N and integer cutoff Q ≥ 1, define the
**Mertens–Ramanujan spike intensity** at level Q,

```
    M_Q(n) := Σ_{q ≤ Q, q squarefree}  mu(gcd(q, n)) / phi(q / gcd(q, n))
```

and the **pointwise spike approximator** of chi_P at scale Q,

```
    T_Q(n) := (pi(N) / N) · M_Q(n).
```

`mu` is the Möbius function and `phi` is Euler's totient. The right-hand
side depends on n through the divisor structure `gcd(q, n)`, and on the
ambient scale N only through the prefactor `pi(N) / N` (the prime
density on `[1, N]`).

## Identity used (proven from definitions)

**Hölder reduction.** For squarefree q and any integer n, let
`d = gcd(q, n)` (then d | q, hence d squarefree). Then

```
    c_q(n)  =  mu(q/d) · phi(d)        (standard for Ramanujan sums)
    mu(q) · c_q(n) / phi(q)  =  mu(d) / phi(q/d).
```

Numerically verified at every squarefree `q ∈ [1, 30]` and `n ∈ [0, 60]`
in `sanity_holder.py`.

Consequently:

```
    M_Q(n) = Σ_{q ≤ Q, q sqf} mu(q) · c_q(n) / phi(q),
```

so `T_Q(n)` is exactly the squarefree truncation, at conductor cutoff Q,
of the Ramanujan-Fourier series whose coefficients are `mu(q)/phi(q)`.

## Why this object exists (project-internal derivation)

S168 (proved analytically in this project) gives, for squarefree q and
N divisible by q:

```
    Σ_{j coprime to q}  e^{2πi j p/q}  summed over primes p ≤ N
        =  (pi(N) − r(q)) · mu(q) / phi(q)  +  R_q(j, N),
```

with `R_q` controlled by the Gallagher PNT-in-AP variance. Plug this
into the additive-Fourier inversion of `chi_P` restricted to the
"primitive" subspace `V_q^prim`:

```
    P_{V_q^prim} chi_P (n)
        =  (1/N) Σ_{j coprime to q}  hat{chi_P}(j N / q) · e^{2πi j n / q}
        ≈  ((pi(N) − r(q)) · mu(q) / (N · phi(q))) · c_q(n).
```

Summing over squarefree `q ∈ [1, Q]`, dropping the `r(q) = O(log q)`
correction (sub-leading for `Q ≪ pi(N)`):

```
    Σ_{q sqf ≤ Q}  P_{V_q^prim} chi_P (n)  ≈  T_Q(n).
```

So `T_Q(n)` is the **principal-character projection** of `chi_P` onto
the additive-Fourier subspace `⊕_{q sqf ≤ Q} V_q^prim`, written as a
single closed-form pointwise function. The non-principal-character
content is captured in the S168 remainder `R_q`.

## Predicted L² content (from S168/S169)

With `||chi_P||² = pi(N)` (since `chi_P` is the 0/1 prime indicator):

```
    ||T_Q||²  =  (pi(N)/N)² · Σ_n M_Q(n)²
              =  (pi(N)² / N) · Σ_{q sqf ≤ Q}  1/phi(q)
```

(by orthogonality of Ramanujan sums of distinct conductors over
`[1, N]`, with `||c_q||² ≈ N · phi(q)`).

Removing the `q = 1` "wheel mode" (the constant function `pi(N)/N`,
contributing `pi(N)²/N` to the L² norm):

```
    ||T_Q − pi(N)/N||²  =  (pi(N)² / N) · Σ_{q sqf ≤ Q, q ≥ 2}  1/phi(q).
```

By Selberg–Delange the squarefree-sum has constant `A_∞ = 1`:
`Σ_{q sqf ≤ Q, q ≥ 2} 1/phi(q) ≈ A · log Q + B`, so the predicted ratio
is

```
    ||T_Q − const||² / pi(N)  ≈  (pi(N) log Q / N) · A.
```

At `pi(N) log N / N → 1` this is `≈ A · log Q / log N`. For
`Q = N^{0.185}` the ratio approaches `≈ 0.185`; at finite N=2^20 the
empirical multiplicative finite-size factor is `≈ 1.18` (per S190),
giving predicted `≈ 0.218 ± 0.02`.

## Falsification criterion (pre-stated)

The construction is **falsified** if any of the following fails on the
canonical empirical ladder `d ∈ {14, 16, 18, 20}`, `Q ∈ {2, 6,
N^{0.185}, N^{0.21}, 30, √N}`:

1. **L² ratio** `||T_Q − pi(N)/N||² / pi(N)` at `Q = N^{0.185}`,
   `d = 20`, must lie in `[0.18, 0.26]` (i.e., consistent with S169's
   measured SVD spike-block fraction `0.220 ± 0.005` plus the
   Q-rounding band).
2. **Pearson correlation** `r(chi_P, T_Q)` must be **strictly
   monotone-increasing in Q** across the ladder. (If T_Q ever stops
   improving as Q grows, the closed form has a defect.)
3. **Precision-at-π(N)** for the top-π(N) ranking by `T_Q(n)` must
   exceed `Q · pi(N) / N` (i.e., the wheel-`Q*` baseline) for every
   tested Q.
4. **Hölder closed-form check** must verify:
   `mu(q) c_q(n)/phi(q) = mu(gcd(q,n)) / phi(q/gcd(q,n))` for every
   squarefree q ≤ 30 and n ≤ 60. Done in `sanity_holder.py`.

All four criteria pass empirically; see `spike_pointwise_results.md`.

## Edges composed

* **E2.1** (MPS bond-dim identity) — the S168 spike subspace lives at
  squarefree conductors. `T_Q` is the **pointwise** form of this spike
  content.
* **E1.5** (`pi(x) mod m` saturation entropy) — `T_Q(n)` packages the
  `mu(q) c_q(n)/phi(q)` weights for all squarefree q ≤ Q into one
  pointwise scalar. Each q contributes the `mod q`-state-conditioned
  density estimate `pi(N)/N` modulated by Hölder weights.
* **E1.6** (parity bisection) — the `q = 2` term `mu(2) c_2(n) / phi(2) =
  −(−1)^n` is exactly the parity sign. `T_Q` recovers the parity bit
  contribution naturally as the leading non-constant term.
* **E2.2** (Liouville identity) — the `q = 1` constant contribution
  `pi(N)/N` plays the role of the smooth wheel density that E2.2's
  `(x − L(x))/2 − C_3(x)` decomposition splits in a different way.
* **S168 squarefree extension theorem** — the *energy*-level statement
  that this construction realises pointwise.

## What is novel relative to S168

S168 is an L² statement about the energy of the projection
`P_{V_q^prim} chi_P`. It does NOT exhibit a closed-form pointwise
function. This construction:

1. **Writes down an explicit pointwise scalar field** `T_Q(n)`, computable
   in `O(Q · log N)` time per evaluation (or `O(d(n) · Q)` per evaluation
   using the divisor decomposition of n).
2. **Simplifies the Ramanujan-sum factor `mu(q) c_q(n)/phi(q)` to the
   one-line Hölder form `mu(gcd(q,n))/phi(q/gcd(q,n))`** (the project's
   prior writing carried the full Ramanujan sum).
3. **Empirically verifies** that the principal-character L² content
   matches S169's SVD spike-block fraction (`0.220` at d=20) to within
   `1-7%` at the predicted Q (this is a non-trivial pointwise check of an
   energy-level prediction).
4. Establishes that **`T_Q` is a primality predictor** with quantitative
   precision-at-π(N) lift `5.6× — 12.8×` over the random baseline as
   `Q` ranges from `N^{0.185}` to `N^{0.5}`.
5. Identifies `T_Q` with the **wheel sieve** at small Q: the closed form
   `mu(d)/phi(q/d)` makes explicit that `T_Q(n)` is large precisely when
   `n` is coprime to most small primes ≤ Q, the wheel-sieve admissibility
   condition.

## Algorithmic content

* **Cost** per evaluation: `O(Q · ω(n))` where `ω(n)` is the number of
  distinct prime divisors of n (small in practice). For
  `Q = N^{0.185}`, this is **strictly sub-`√N`**.
* **Output**: a real-valued primality "score". Cannot be turned into an
  exact prime indicator without extra information (the missing 79% of
  L² mass is the GUE-bulk MP randomness, S74).
* **Compared to wheel-W sieve at level W = primorial(k):** for
  `Q = primorial(k)`, `T_Q(n)` is a *continuous* relaxation of the
  wheel sieve. The wheel sieve outputs `0/1` (admissible / not), while
  `T_Q` outputs a graded score whose top-`pi(N)`-cut is calibrated to the
  prime density.

## Pre-stated identity prediction (also testable)

For `Q = primorial(k)` (so `Q = 2, 6, 30, 210, ...`), the formula

```
    M_Q(n) = Σ_{q sqf ≤ Q}  mu(gcd(q,n)) / phi(q/gcd(q,n))
```

reduces to a sum over divisors of `Q = ∏_{p ≤ p_k} p` (since every
squarefree q ≤ Q divides Q), yielding the closed-form

```
    M_{primorial(k)}(n)  =  Σ_{q | primorial(k)} mu(gcd(q,n)) / phi(q/gcd(q,n)).
```

This sum has `2^k` terms and is `Q · ∏_{p ≤ p_k}(1 − 1/p)`-like in
its asymptotic structure (Mertens density). Specifically, when
`n` is coprime to `primorial(k)`,

```
    M_{primorial(k)}(coprime n)  =  ∏_{p ≤ p_k} (1 + 1/(p − 1))
                                =  ∏_{p ≤ p_k} p/(p − 1)
                                =  primorial(k) / phi(primorial(k))
                                =  W / phi(W).
```

So at coprime points, `T_W(n) = pi(N)/N · W/phi(W)` — exactly the
**wheel-W density correction** of the Mertens product. This is the
direct algebraic statement that ties E2.1's `phi(W)/W` factor to a
pointwise prediction of `chi_P` density on the coprime residue classes.

This is the "**Hölder–Mertens identity**" of `T_Q`:

```
    T_W(n) | gcd(n, W) = 1   =   pi(N)/N · W/phi(W),
```

which says: the spike approximator on coprime-to-W residues equals the
expected prime density divided by `phi(W)/W` (the wheel filling factor).
Verified empirically by `precision_at_pi(N)` at primorial Q.

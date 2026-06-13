# Möbius bisection of π(x) — definition

## Setting

Let `λ(n) = (-1)^Ω(n)` (Liouville), `μ(n) = (-1)^ω(n)·1[n squarefree]`
(Möbius), `1_sqf(n) = μ(n)²`. Define cumulative summatories:

* `L(x) = Σ_{n ≤ x} λ(n)`           (Liouville summatory)
* `M(x) = Σ_{n ≤ x} μ(n)`           (Mertens function)
* `Q(x) = Σ_{n ≤ x} 1_sqf(n)`        (squarefree count)
* `π(x)` = prime-counting function

## The 4-cell decomposition

Partition `[1, x]` by the pair `(squarefree-ness, Ω-parity)`:

|                 | Ω even   | Ω odd    | row sum |
|-----------------|----------|----------|---------|
| **squarefree**  | `S_e(x)` | `S_o(x)` | `Q(x)`  |
| **non-sqfree**  |`NS_e(x)` |`NS_o(x)` |`x − Q(x)`|

For a squarefree `n`, `Ω(n) = ω(n)` and `λ(n) = μ(n)`. Hence `M = S_e − S_o`
and `Q = S_e + S_o`. For a non-squarefree `n`, `μ(n) = 0`. So
`L − M = NS_e − NS_o` (because `λ(n) = μ(n)` on squarefree, `μ = 0` off it).
Solving the four linear equations:

```
   S_e  = (Q + M) / 2          NS_e = ((x − Q) + (L − M)) / 2
   S_o  = (Q − M) / 2          NS_o = ((x − Q) − (L − M)) / 2
                                    = (x − Q − L + M) / 2
```

## The two bisections

The squarefree-Ω-odd cell `S_o(x)` decomposes as
`S_o(x) = π(x) + C_3*(x)`, where

```
   C_3*(x) := #{n ≤ x : sqfree, composite, ω(n) odd, ω(n) ≥ 3}.
```

(Squarefree integers with ω = 1 are exactly the primes; ω = 3, 5, 7, ...
are squarefree composites with odd ω.) Therefore:

> **Möbius bisection**
> `π(x) = (Q(x) − M(x)) / 2 − C_3*(x)`.

The Liouville bisection (E2.2) decomposes the larger Ω-odd cell
`N_o(x) = S_o(x) + NS_o(x)`:

```
   π(x) = (x − L(x)) / 2 − C_3(x),
   C_3(x) := #{n ≤ x : composite, Ω(n) odd, Ω(n) ≥ 3}
           = C_3*(x) + NS_o(x).
```

Subtracting the two bisections and using `S_o = π + C_3*` and
`N_o = π + C_3` gives the bridge identity:

> **Bridge identity (NS_o closed form)**
> `NS_o(x) = (x − Q(x) − L(x) + M(x)) / 2`.

Equivalently `C_3(x) − C_3*(x) = NS_o(x)`.

The bridge identity admits a one-line analytic proof. The indicator of
`{not sqfree, Ω odd}` is `(1 − μ²(n)) · (1 − λ(n)) / 2`, so summing,

```
   NS_o(x) = (1/2) Σ_{n ≤ x} (1 − μ²(n))(1 − λ(n))
           = (1/2) (x − Q(x) − L(x) + Σ μ²(n)λ(n)).
```

For squarefree `n`, `μ²(n) = 1` and `λ(n) = μ(n)`; for non-squarefree `n`,
`μ²(n) = 0`. Hence `Σ_{n ≤ x} μ²(n)λ(n) = Σ_{n ≤ x} μ(n) = M(x)`, giving
`NS_o = (x − Q − L + M) / 2`. ∎

## Edges composed

* **E2.2** Liouville bisection `π = (x − L)/2 − C_3`.
* **E1.6 / E2.10** parity bisection `π(x) mod 2 = A(x) mod 2 ⊕ C_3(x) mod 2`,
  with `L(x) mod 2 = x mod 2` (the parity trap).
* The Möbius `M = S_e − S_o` decomposition is folklore (Mertens 1874);
  this is its π(x)-bisection consequence — not previously recorded as a
  project edge.

## Falsification criteria (pre-stated)

For all integers x ≤ N (set N = 10⁶ for production):

* **F1** `π(x) = (Q(x) − M(x))/2 − C_3*(x)` integer-exact.
* **F2** `NS_o(x) = (x − Q − L + M)/2` integer-exact.
* **F3** `C_3(x) − C_3*(x) = NS_o(x)` integer-exact.
* **F4** Parity: `((x − L)/2 ⊕ (Q − M)/2) ≡ NS_o(x) (mod 2)`.
* **F5** Möbius parity bisection: `π(x) ≡ (Q − M)/2 + C_3*(x) (mod 2)`.
* **F6** Liouville parity bisection (sanity): `π(x) ≡ (x − L)/2 + C_3(x) (mod 2)`.

If any falsifier fires, the construction is wrong.

## What this construction is and is not

* It is a **textbook-level identity** rearranging Möbius / squarefree /
  Liouville summatories. A working analytic number theorist would derive
  it in well under an hour.
* It is **not previously recorded in the project** — neither the Möbius
  bisection itself, nor the bridge identity `NS_o = (x − Q − L + M)/2`,
  nor the 4-cell decomposition, appear in `EDGES.md` or `CLOSED_PATHS.md`.
  E2.2 records only the Liouville-side bisection; E1.6 records only the
  Liouville-side parity decomposition.
* It is **not an algorithmic improvement**. `Q(x)` is `Õ(√x)`, but
  `M(x)` and `C_3*(x)` are both `O(x^{2/3})` to compute by best known
  methods, and `(Q − M)/2 ~ 0.304 · x` is dominated by `C_3* ~ 0.304 · x`
  with `π ~ x/log x` as the small "needle" — the same C-circular obstacle
  as E2.2's `C_3`.
* It **is** a structural unification: π(x) admits at least two parity
  bisections (Liouville-side over [1, x]; Möbius-side over the
  squarefree subset), and the two are linked by the explicit closed
  form for the count of non-squarefree Ω-odd integers.

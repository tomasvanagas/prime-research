# Spike-fraction specificity: is S168's 21% a primality fact or a generic Fourier feature?

**Construction date:** 2026-04-28 (Session 200, paradigm-shift mode).
**No cross-domain technique imported.**

## Object

Define the **squarefree-Q spike fraction** of an arithmetic indicator
`f : [1, N] → {0, 1}`:

```
   Spike(f, Q, N)  :=  Σ_{q sqf, 2 ≤ q ≤ Q}  ||P_{V_q^prim} f||² / N
                       ────────────────────────────────────────────────
                                          ||f||² / N
```

where `V_q^prim` is the primitive-additive-character subspace at
conductor `q`, i.e. the span of `e^{2πi a n / q}` for `a` coprime to `q`.

Computationally, with `F_q[r] = Σ_{n ≤ N, n ≡ r (mod q)} f(n)`:

```
   ||P_{V_q^prim} f||²  =  (1/N) Σ_{a coprime to q}  | Σ_r F_q[r] e^{-2πi a r / q} |²
```

## Edges composed

* **E2.1** (S168 squarefree-q formula): `E(q, N; χ_P) ≈ μ(q)²·(π(N))² / (φ(q) N)`.
  S168/S190 conclude that `Σ_{q sqf ≤ N^{0.185}} E(q, N; χ_P) ≈ 0.21 · π(N)`.
* **E2.2** (Liouville parity identity): `π(x) = (x − L(x))/2 − C₃(x)`, where
  `(x − L(x))/2 = #{n ≤ x : Ω(n) odd}` and `C₃(x) = #{n ≤ x : Ω(n) odd, ≥ 3}`.
* **E2.10** (free identity `L(x) mod 2 = x mod 2`): the Liouville-parity
  decomposition splits π into two MORE pseudorandom components.

## Question

S168/S190 establish the 21% spike fraction for `χ_P`. Is this fraction:

  - **(P)** a specific consequence of "primes have no support on residues
    divisible by `q`" (the principal-character / coprime-support
    mechanism), OR
  - **(G)** a generic feature of arithmetic indicators with density
    `~ 1 / log N`?

## Hypothesis (pre-stated)

**(P) is correct.** Specifically:

- `Spike(χ_P, Q, N) ≈ 0.21` at `Q = N^{0.185}` (control, S168/S190).
- `Spike(χ_{Ω odd}, Q, N) → 0` as N → ∞ (Liouville indicator: full
  residue support, no principal-character mechanism).
- `Spike(χ_{Ω = 2}, Q, N)` (semiprime indicator): small (O(1/log N)).
- `Spike(χ_{Ω = 3}, Q, N)`: small.
- `Spike(μ², Q, N)` (squarefree indicator): non-zero but distinct
  fingerprint (residue bias is `(1 − 1/p²)`, not `χ_P`'s `1/(p − 1)`).

## Falsification criteria (committed before run)

| Criterion | If TRUE, hypothesis (P) | If FALSE, hypothesis (G) |
|---|---|---|
| `Spike(χ_P) > 5 · Spike(χ_{Ω odd})` at d=16, Q=8 | survives | refuted |
| `Spike(χ_{Ω odd}) < 0.05` at d=16, Q=8 | survives | refuted |
| `Spike(χ_P) > 3 · Spike(χ_{Ω=2})` at d=16, Q=8 | survives | refuted |
| `Spike(χ_P)` decreases from d=14 → d=16 (finite-d gap) | survives | refuted |

If all four hold: hypothesis (P) confirmed — **the spike fraction is a
quantitative arithmetic signature of "no small prime factors" support
structure, not a pseudorandomness feature**.

If the contrast `χ_P / χ_{Ω odd}` ratio is < 2: hypothesis (P) rejected
— the spike fraction is generic and S168's 21% has been over-interpreted
as primality-specific.

## Why this matters

The 21% spike fraction is currently anchored at S190 as one half of a
two-piece decomposition (`π(N) = 0.21 · π(N) [PNT-in-AP at sqf q ≤
N^{0.185}] + 0.79 · π(N) [GUE-bulk MP]`). This decomposition is treated
as a structural fact about `χ_P`. If hypothesis (G) holds, the decomposition
is a generic Fourier identity unrelated to primality, and the structural
weight of the S168/S190 result drops substantially. If (P) holds, the
decomposition is a quantitative refinement of the fact that primes are
the sole arithmetic indicator with the "(π(N))²/(φ(q) N)" principal-
character contribution, which IS primality-specific.

## Algorithmic implication (limited but real)

If (P) holds: the T_Q-style polylog spike approximator (S191) does NOT
extend to `χ_{Ω odd}` (Liouville indicator) or to `χ_{Ω = 2}` (semiprime
indicator). The polylog primality-prefilter idea is anchored in the
specific arithmetic structure of "coprime support", and any future
"Liouville-prefilter" / "semiprime-prefilter" of the same flavour is
ruled out by this construction.

If (G) holds: similar polylog prefilters DO exist for chi_{Ω odd} and
chi_{Ω = 2}. This would be a substantial expansion of the S191
algorithm's scope.

## Files

* `spike_indicator_specificity.py` — main script.
* `spike_indicator_specificity_results.md` — empirical outcome.
* `spike_indicator_specificity.json` — raw per-(f, Q, d) numerics.
* `definition.md` — this file.

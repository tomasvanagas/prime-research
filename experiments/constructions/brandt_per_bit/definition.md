# C3 — Brandt obstructions × per-bit difficulty: formal definitions

## Composition

This construction composes two project edges:

- **E5.8** [Brandt-class barrier]. Brandt 2024 (TCC, ePrint 2024/687) proves
  `MKtP not in DTIME[O(n)]` by a TRAVERSE diagonalisation that requires a
  Kt-random hard string of growing length. Four obstructions (O1-O4)
  prevent the technique from extending to fixed natural functions like
  `π(x) mod 2`:
  (O1) the hard string is oracle-dependent, not fixed;
  (O2) the contradiction uses self-referential Kt on both sides;
  (O3) the only ingredient bypassing the black-box barrier is the
       1-Kt-randomness of Chaitin Ω prefixes (no analog for fixed total
       Boolean functions);
  (O4) the bound is on uniform time, not on circuits.
- **E1.3** [4-bit difficulty boundary]. The per-bit Fourier / influence
  / R-correlation profile of `π(x)` shows a sharp transition near bit
  position J ≈ 0.5N: bits J in the **easy zone** `[0.6N, N-1]` (top of
  the integer) are highly correlated with the smooth `R^{-1}(n)`
  approximation; bits J in the **hard zone** `[0, 0.4N]` (bottom) are
  essentially random under the project's pseudorandomness battery.

## The objects

### Per-bit prime-counting function

For each integer J ≥ 0 and each N ≥ 1 we define

```
    π_J : N → {0,1}
    π_J(x) := bit J of π(x)  (binary representation, J=0 is LSB)
```

The associated finite-N truth-table is

```
    s_J^{(N)} := ( π_J(0), π_J(1), ..., π_J(2^N - 1) ) ∈ {0,1}^{2^N}.
```

Note `π_0` = `π(x) mod 2` is exactly Brandt's MKtP-encoded target from
session 51 (`brandt_mktp/`). For J ≥ 1, `π_J` is a NEW Boolean function,
not previously tested under the bounded-Kt framework.

### The per-bit family

The composition object is the **family**

```
    F_N := { s_J^{(N)} : J = 0, 1, ..., N-1 }.
```

For fixed N, this is a finite collection of N truth-tables, each of
length 2^N. Notes:

- For sufficiently large J, `π(x) < 2^J` for all `x < 2^N`, so
  `π_J ≡ 0` and `s_J^{(N)} = 0^{2^N}`. The first such J is
  `J*(N) := ⌈log₂(π(2^N - 1) + 1)⌉`. For N ≤ 7 this is roughly
  `N - 2`.
- For small J, `s_J^{(N)}` looks pseudorandom (E1.9 measures); the
  empirical question this construction asks is whether the bounded-Kt
  simulator from `brandt_mktp/` distinguishes the easy / hard halves.

## The relationship to π(x)

Knowing `s_0^{(N)}, s_1^{(N)}, ..., s_{N-1}^{(N)}` simultaneously
recovers the full `π(x)` table. So the per-bit family is an
**informationally complete** decomposition of `π` on `[0, 2^N)`.

The Brandt-extension question then sharpens to:

> Is there a J such that `π_J` is bounded-Kt-non-random (compressible by
> our 3-bit-per-op stack VM at `L_MAX=12`) while `π_0` is not?
>
> If yes: a Brandt-style TRAVERSE on the per-bit family could in
> principle exploit J-stratified Kt-randomness — the "minimal weakening
> of O1" the C3 spec asks about.
>
> If no: per-bit decomposition does not provide additional structure
> beyond the wholesale bound on `π_0` (DUPLICATE-PLUS of E5.8).

## Bounded-Kt simulator (imported from `brandt_mktp/`)

We re-use the 3-bit-per-instruction stack VM with operations
`{PUSH 0, PUSH 1, DUP, FLIP, RPT2, RPT4, ALT, HALT}`, `L_MAX = 12`,
`T_MAX = 4096`. For a target string `s` we compute

```
    Kt_b(s) := min { L + ⌈log₂ t⌉  :  Π is L-bit, U(Π) outputs s in t ≤ T_MAX steps }
```

over `L ∈ {3, 6, 9, 12}`. Saturation `Kt_b = INF = 61` indicates no program
in the VM produces `s`. Compression `Kt_b ≪ |s|` indicates structural
regularity expressible in the VM's primitives.

## Falsification (pre-stated)

The construction tests three pre-stated outcomes; exactly one of these
will hold after running the experiment:

- **F1 (collapse / DUPLICATE-PLUS).** `Kt_b(s_J^{(N)}) = 61`
  (saturated) for **every** `J ∈ [0, N-1]` and every tested `N`. Reading:
  the bounded VM cannot distinguish bit positions; per-bit decomposition
  gives no new structure under bounded-Kt; C3 reduces to E5.8 wholesale.
  (Closure: mode E, DUPLICATE-PLUS of S51.)

- **F2 (J-stratified, matching E1.3).** For each tested N:
  - All `J < 0.5·N` (the *hard zone*): `Kt_b(s_J^{(N)}) = 61` (saturated);
  - At least one `J ≥ 0.5·N` (the *easy zone*): `Kt_b(s_J^{(N)}) < 61`
    (compressible by the VM).
  Reading: the bounded VM **does** distinguish the easy / hard halves.
  This sharpens E1.3 from a Fourier-spectrum statement to a Kt-complexity
  statement and gives a per-bit family with non-trivial structure under
  Brandt's specific complexity measure. (Composition genuinely novel; the
  per-bit-Kt-stratified family did not exist in the project before.)

- **F3 (anomalous transition).** F1 and F2 both fail: the J-Kt-saturation
  pattern is non-monotonic, or the easy/hard cut differs from E1.3's
  0.5N boundary by more than 1 bit position. Reading: bounded-Kt sees
  structure that R-correlation / Fourier do **not**, or vice versa —
  potentially a NEW orthogonal pseudorandomness-vs-compressibility edge.

## Edge IDs cited

- E5.8 (Brandt-class barrier)
- E1.3 (per-bit difficulty profile)
- E1.9 (pseudorandomness of π(x) mod 2 under 35+ measures)
- E5.3 (PRIMES TC⁰ open) — the problem the composition lives near

## Files in this directory

- `brandt_per_bit.py` — bounded-Kt evaluator + per-bit table builder
  + outcome verdict (F1 / F2 / F3).
- `definition.md` — this file.
- `brandt_per_bit_results.md` — pre-stated falsifiers + run output +
  verdict + structural analysis of which obstruction (O1-O4) shifts.

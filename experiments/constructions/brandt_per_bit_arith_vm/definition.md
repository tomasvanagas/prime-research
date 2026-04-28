# C3.a — Arithmetic-primitive bounded-Kt VM — definition

## Composition

This construction is the explicit successor of C3 (S105) proposed in
`brandt_per_bit_results.md`. It composes the same two edges:

- **E5.8** — Brandt 2024 is structurally welded to MKtP. The structural
  obstructions O1 / O2 / O3 / O4 to extending Brandt's TRAVERSE
  diagonalisation to fixed natural Boolean functions.
- **E1.3** — Per-bit difficulty of `p(n)` (and analogously of `π(x)`)
  has a sharp ~4-bit sigmoid boundary. For `π(x)` over `x ∈ [0, 2^N)`,
  the bits of `π(x)` partition into:
    * **hard** (LSB side, `J < ⌈N/2⌉`) — pseudorandom parities;
    * **easy** (MSB side, `J ≥ ⌈N/2⌉`) — match `Li(x)` smoothly;
    * **trivial** (`J ≥ J*(N)` where `J*(N) := ⌈log₂(π(2^N) + 1)⌉`) —
      identically zero because `π(2^N) < 2^J`.

S105 found that the original 3-bit-per-op stack VM resolves only the
trivial cut at `J*(N) ≈ N − log₂ N`, **not** E1.3's smooth/oscillatory
cut at `0.5N`. The bounded-Kt simulator was structurally blind to the
arithmetic content distinguishing easy bits from hard bits.

C3.a tests whether **augmenting the VM with arithmetic primitives**
(LOG2, LI_APPROX, DIV_LOG, GEO_SUM — the spec items proposed in S105's
successor block) shifts the bounded-Kt cut from `J*(N)` toward
`⌈N/2⌉`. The hypothesis is: if the VM can compute `Li(x)`, it should
compress the easy-zone bits `J ∈ [⌈N/2⌉, J*(N))` because those are
exactly the bits where `bit_J(π(x)) ≈ bit_J(Li(x))` to within
finite-x noise.

## Object: extended bounded-Kt simulator

A 4-bit-per-op stack VM with **integer stack** (separate from the
emitted bit stream) and the following 16 ops:

| Bits   | Mnemonic   | Effect                                                                  |
|--------|------------|-------------------------------------------------------------------------|
| 0x0    | PUSH0      | Push 0 onto integer stack S.                                            |
| 0x1    | PUSH1      | Push 1.                                                                 |
| 0x2    | PUSH_N     | Push current emit-count `len(out)` (the "step counter" / x-input).      |
| 0x3    | INC        | TOS += 1 (saturating at INT_CAP).                                       |
| 0x4    | ADD        | Pop a, b; push min(a + b, INT_CAP).                                     |
| 0x5    | SUB        | Pop a, b; push max(b − a, 0).                                           |
| 0x6    | MUL        | Pop a, b; push min(a · b, INT_CAP).                                     |
| 0x7    | SHR1       | TOS >>= 1.                                                              |
| 0x8    | LOG2       | TOS ← max(TOS, 1).bit_length() − 1   (i.e. ⌊log₂ TOS⌋, ⌊log₂ 0⌋ := 0).  |
| 0x9    | LI_APPROX  | Pop a; push ⌊a / log a⌋ if a ≥ 2 else 0.   (R^{-1} kernel.)             |
| 0xa    | DIV_LOG    | Pop a, b; push ⌊b / log max(a, 2)⌋.                                    |
| 0xb    | GEO_SUM    | Pop a; push 1 + a + a² + ... terminated when partial sum > T_MAX.       |
| 0xc    | DUP        | Duplicate TOS.                                                           |
| 0xd    | EMIT_S     | Pop a; emit '1' if a & 1 else '0'.                                      |
| 0xe    | EMIT0      | Emit '0' (no stack change).                                              |
| 0xf    | HALT       | Stop.                                                                    |

Programs of length L (multiple of 4) loop: `pc` wraps to 0 when it
reaches L. Caps: `INT_CAP = 10⁹`; `T_MAX = 4096`; `L_MAX_EXT = 24`
(6 instructions). Bounded Kt:
```
   Kt_b'(t) := min { L + ⌈log₂ steps⌉ : prog of length L emits t in ≤ T_MAX steps }
```
or `INF = 4·L_MAX_EXT + ⌈log₂ T_MAX⌉ = 109` if no program produces `t`.

## Per-bit family s_J^(N)

Same as C3:
```
   s_J^(N) := ( bit_J(π(0)), bit_J(π(1)), ..., bit_J(π(2^N − 1)) ) ∈ {0,1}^{2^N}
```
with bit 0 = LSB. We compute `Kt_b'(s_J^(N))` for `N ∈ {3, 4, 5, 6}`,
`J = 0 .. N−1`. All cells from C3-S105 are re-run; we directly compare
`Kt_b` (original VM) to `Kt_b'` (extended VM) on the SAME bit strings.

## Pre-stated falsifiers

Pick exactly one of the three based on the empirical Kt_b' cut location:

- **F1 — full shift to E1.3 boundary.** For every tested N ≥ 4, the
  extended-VM Kt_b' is < INF for **every** `J ≥ ⌈N/2⌉` (the entire
  E1.3 easy zone), and Kt_b' = INF for every `J < ⌈N/2⌉`. The cut
  location moves from `J*(N) ≈ N − log₂ N` to `⌈N/2⌉` exactly.
  Implication: the smooth/oscillatory transition E1.3 detects via
  Fourier weight becomes bounded-Kt-visible once the VM has R^{-1}-
  kernel primitives. Genuine new edge content.

- **F2 — no shift.** For every N ≥ 4, Kt_b' = INF for every
  `J ∈ [⌈N/2⌉, J*(N))` (i.e. the easy-zone bits between E1.3's
  prediction and the trivially-zero floor remain saturated). The
  arithmetic primitives don't help compress finite-x π bits because
  the π − Li residual is incompressible at this resolution.
  Implication: refines E1.3 — the smooth/oscillatory split is a
  *Fourier-spectrum* fact, not a *bounded-Kt-program* fact, even with
  R^{-1}-aware primitives. DUPLICATE-PLUS of E5.8 / C3.

- **F3 — intermediate hierarchy.** For at least one tested N ≥ 4,
  Kt_b' < INF for some `J ∈ [⌈N/2⌉, J*(N))` but Kt_b' = INF for
  another `J' ∈ [⌈N/2⌉, J*(N))`. The cut sits **strictly between**
  `⌈N/2⌉` and `J*(N)`. Implication: a new VM-richness-indexed
  hierarchy of cut locations. Identify the threshold cell(s) and what
  programs realize them.

## Why this is a valid composition (not a duplicate)

C3 / E5.8 closes Brandt's structural argument for the per-bit family
at the **logical** level (O1..O4 still apply). C3.a is **orthogonal**:
it asks an empirical / quantitative question about a different object —
the bounded-Kt complexity of the same family under a different VM —
and the structural argument (which doesn't depend on the VM) is
unchanged regardless of the empirical outcome. The CONSTRUCTION is the
extended VM and the empirical Kt_b' table; the COMPOSITION mechanism
is "use E1.3's bit-structure to choose targets, use E5.8's bounded-Kt
machinery to evaluate them, use arithmetic primitives to bridge".

## Edges cited

- **E5.8** — Brandt 4-obstruction argument (inherited from C3; the
  structural conclusion is unchanged by VM choice).
- **E1.3** — Per-bit difficulty profile (the targets `s_J^(N)` come
  from this edge; the cut-location question is a quantitative
  refinement).
- **E1.5** — Information rate `h_2(π(X)/X)` per step (the "easy zone"
  exists exactly because the smooth predictor `Li⁻¹` saturates the
  per-step information rate; no contradiction with E1.5).
- **E5.3** — PRIMES TC⁰ (the open problem the composition lives near
  — neither F1 nor F2 nor F3 produces a TC⁰ circuit; the construction
  inherits the structural distance from E5.3).

## Falsifier outcomes — what each produces

| Outcome | Status              | Filing                                             |
|---------|---------------------|----------------------------------------------------|
| F1      | Genuine new edge    | `novel/arith_vm_pi_bit_compression.md` + EDGES row |
| F2      | DUPLICATE-PLUS      | CLOSED_PATHS row, refine E1.3 inline               |
| F3      | New empirical fact  | EDGES.md inline refinement of E1.3                 |

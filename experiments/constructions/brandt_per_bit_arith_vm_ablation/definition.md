# C3.a.iv — Arithmetic-primitive ablation of the C3.a bounded-Kt VM

## Composition

C3.a.iv extends C3.a (S150) by ablating individual arithmetic primitives
to identify the **causally responsible** primitive(s) for the bounded-Kt
cut shift documented in S150. It composes:

- **E5.8** — Brandt 4-obstruction structural barrier (inherited from C3 / S105).
- **E1.3** — Per-bit difficulty cut at `J = ⌈N/2⌉` (the per-bit family is the target).
- **C3.a / S150** — Extended VM with arithmetic primitives `{LOG2, LI,
  DIV_LOG, GEO_SUM}`. S150 found `bit_J(π)` cells in the easy zone
  compress at `L_max ∈ {24, 28}`, with optimum programs that variously
  use `LI` (N=4 J=2, N=5 J=3, N=6 J=4), `DIV_LOG` (N=3 J=1), and
  `GEO_SUM` (N=3 J=0).

The composition mechanism: re-run S150's exact target battery under the
same C inner loop with one or more arithmetic primitives **forbidden**.
Compare per-cell `Kt_b'` shifts to the S150 baseline. The empirical
question is whether one primitive (LI, the literal `Li(x) = x / log x`
kernel) carries the entire easy-zone shift, or whether the four
primitives are jointly required at different `(N, J)` cells.

## Object

A bit-mask-parametrised version of S150's batch enumerator. The 16-bit
`forbidden_ops_mask` selects ops to skip; any program containing a
forbidden op is rejected pre-simulation. All other VM mechanics are
identical to S150 (T_MAX = 4096, INT_CAP = 10⁹, 4-bit ops, looping pc).

Six ablation conditions:

| name           | forbidden ops              | tests                                           |
|----------------|----------------------------|-------------------------------------------------|
| `baseline`     | (none)                     | sanity vs S150                                  |
| `drop_LOG2`    | LOG2                       | LOG2 necessary?                                 |
| `drop_LI`      | LI                         | LI necessary? (predicted: yes for easy zone)    |
| `drop_DIV_LOG` | DIV_LOG                    | DIV_LOG necessary?                              |
| `drop_GEO_SUM` | GEO_SUM                    | GEO_SUM necessary?                              |
| `only_LI`      | LOG2, DIV_LOG, GEO_SUM     | LI sufficient? (predicted: yes for easy zone)   |

Each ablation reproduces the per-bit `Kt_b'` table for `N ∈ {3, 4, 5, 6}`
and `J = 0..N−1`. The full ablation grid takes ≲ 1 minute at L_max=28.

## Pre-stated falsifiers

S150 baseline shows easy-zone compression at exactly two (N ≥ 4) cells
under L_max=24 and one more under L_max=28:

- N=4 J=2: `ADD, LI, LI, EMIT_S, PUSH_N, PUSH_N` (uses LI)
- N=5 J=3: `EMIT_S, LI, DUP, INC, PUSH_N, PUSH0` (uses LI)
- N=6 J=4: `EMIT_S, PUSH_N, LI, LI, LI, DUP, EMIT_S` (uses LI; L_max=28 only)

We test five mutually-exclusive verdicts against the easy-zone subset
`{(N, J) : ⌈N/2⌉ ≤ J < J*(N), N ≥ 4}`:

- **F1** — `drop_LI` saturates **every** easy-zone cell that baseline
  compresses; `drop_LOG2`, `drop_DIV_LOG`, `drop_GEO_SUM` each preserve
  every easy-zone compression; `only_LI` preserves every easy-zone
  compression. Verdict: **LI is solely necessary and sufficient** for
  the easy-zone shift.

- **F2** — Multiple primitives are jointly necessary. `drop_LI` breaks
  some easy-zone cell *and* at least one of `drop_LOG2`, `drop_DIV_LOG`,
  `drop_GEO_SUM` also breaks at least one easy-zone cell.

- **F3** — `only_LI` leaves at least one easy-zone cell saturated that
  baseline compressed. Verdict: **LI alone does not suffice**; some
  composition with another primitive is required.

- **F4** — Every single-drop preserves every baseline easy-zone
  compression. Verdict: **no single primitive is strictly necessary**;
  redundancy in the VM allows alternative programs.

- **F5 (mixed)** — None of the above clean cases hold; e.g. `drop_LI`
  preserves but `only_LI` breaks (incoherent given monotonicity?).
  Captured in the runner's "mixed" branch.

## Hard-zone informational outcome

S150 compressed two **hard-zone** cells (N=3 J=0 with GEO_SUM; N=3 J=1
with DIV_LOG) at small N. These are below the trivial cut (J* = 3 at
N=3 means J=0, 1 are non-trivial-zero) and the compression is
empirical, not structural. We record per-cell ablation impact for
documentation but the verdict above is anchored on the **easy-zone
N ≥ 4** cells where E1.3's 0.5N-boundary lives.

## Edges cited

- **E5.8** — Brandt structural obstructions, unchanged by ablation.
- **E1.3** — Per-bit difficulty cut, the target family.
- **E1.5** — Information rate; the easy zone is exactly where Li⁻¹
  saturates the per-step rate, providing the algorithmic motivation
  for the LI primitive.
- **E5.3** — PRIMES-in-TC⁰ open problem; the construction is in this
  problem's neighbourhood but does not produce a TC⁰ circuit.

## Falsifier outcomes — what each produces

| Outcome | Status              | Filing                                              |
|---------|---------------------|-----------------------------------------------------|
| F1      | E1.3 inline refine  | Causal-LI fact; CLOSED_PATHS row, EDGES annotation  |
| F2      | E1.3 inline refine  | Joint-primitive necessity fact; CLOSED_PATHS row    |
| F3      | A-grade candidate   | LI insufficient → composition mechanism is novel    |
| F4      | C-grade duplicate   | Redundancy verified; CLOSED_PATHS row               |
| F5      | partial             | Document residual structure                         |

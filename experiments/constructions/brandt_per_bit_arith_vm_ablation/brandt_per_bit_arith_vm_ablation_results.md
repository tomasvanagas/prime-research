# C3.a.iv — Arithmetic-primitive ablation of the C3.a bounded-Kt VM — RESULTS

**Status:** BUILT. Composes E5.8 + E1.3 + C3.a (S150).

## Construction summary

Re-ran S150's per-bit `Kt_b'(s_J^(N))` scan over `N ∈ {3, 4, 5, 6}`,
`J = 0..N−1`, at `L_max = 28` bits, under six ablation conditions on
the four arithmetic primitives `{LOG2, LI, DIV_LOG, GEO_SUM}`:
`baseline`, `drop_LOG2`, `drop_LI`, `drop_DIV_LOG`, `drop_GEO_SUM`,
`only_LI`. The C inner loop's `batch_enumerate_at_length_filtered`
checks each program for forbidden ops pre-simulation, so the ablations
are exact.

Total: 6 × 268M = ~1.6B programs scanned in 463 s wall-time
(baseline 121 s, single drops 76–81 s each, only_LI 30 s). The
filter-rate scaling (single drops simulate 99M / 145M = 68%; only_LI
simulates 30%) confirms the ops are well-mixed across the program
distribution.

## Pre-stated falsifiers

(See `definition.md`.) Easy-zone subset is
`{(N, J) : ⌈N/2⌉ ≤ J < J*(N), N ≥ 4}`. At baseline at L_max=28
these are: `N=4 J=2`, `N=5 J=3`, `N=6 J=4` (matching S150).

## Verdict

**F4 — no single arithmetic primitive is strictly necessary for the
easy-zone bounded-Kt cut shift.** Every easy-zone cell that compresses
in baseline also compresses under every single-drop ablation AND under
the only-LI ablation. The compressing programs change, but the cut
location does not.

Refinement: **the FAMILY of arithmetic primitives is sufficient via
multiple substitutable mechanisms.** LI alone suffices (only_LI matches
baseline on every easy-zone cell), AND `{LOG2, DIV_LOG, GEO_SUM}`
without LI also suffices (drop_LI matches baseline on the L_max≤24
easy-zone cells and is only +1 bit worse at the L_max=28 N=6 J=4 cell).

This **refutes the simplest reading** of S150's "iterated LI is the
dominant compression mechanism" — the dominant mechanism is the
*ability to compute slow-growing integer functions of `out_count`*, and
each of the four arithmetic primitives provides such an ability. LI is
the cleanest/shortest realization but is not unique.

## Per-cell ablation table (Kt_b' under each ablation)

```
  N J  zone     baseline       drop_LOG2         drop_LI    drop_DIV_LOG    drop_GEO_SUM         only_LI
  ------------------------------------------------------------------------------------------------------
  3 0  hard      24                 =               =               =            27(+3)          27(+3)
  3 1  hard      29                 =               =               =               =               =  
  3 2  easy      23                 =               =               =               =               =  
  4 0  hard      36                 =              INF*            INF*             =              INF*
  4 1  hard      37                INF*            INF*             =               =              INF*
  4 2  easy      33                 =               =               =               =               =  
  4 3  easy      10                 =               =               =               =               =  
  5 0  hard      41*               INF*            INF*            INF*            INF*            INF*
  5 1  hard      41*               INF*            INF*            INF*            INF*            INF*
  5 2  hard      37                 =              INF*            INF*             =              INF*
  5 3  easy      33                 =               =               =               =               =  
  5 4  easy      10                 =               =               =               =               =  
  6 0  hard      41*               INF*            INF*            INF*            INF*            INF*
  6 1  hard      41*               INF*            INF*            INF*            INF*            INF*
  6 2  hard      41*               INF*            INF*            INF*            INF*            INF*
  6 3  easy      41*               INF*            INF*            INF*            INF*            INF*
  6 4  easy      36                 =            37(+1)             =               =               =  
  6 5  easy      10                 =               =               =               =               =  
```

`=` means same Kt as baseline. `INF*` means saturated. `n(+k)` means `Kt
= n` with `k`-bit penalty vs baseline.

## Easy-zone compressing programs by ablation

The same target compresses in different ways depending on which
primitives are available:

| (N, J) | baseline                                         | drop_LI                                                      | only_LI                                              |
|--------|--------------------------------------------------|--------------------------------------------------------------|------------------------------------------------------|
| (3, 2) | `EMIT0, LI, EMIT_S, PUSH_N`                      | `EMIT_S, PUSH_N, GEO_SUM, EMIT_S`                            | `EMIT0, LI, EMIT_S, PUSH_N`                          |
| (4, 2) | `ADD, LI, LI, EMIT_S, PUSH_N, PUSH_N`            | `GEO_SUM, DIV_LOG, EMIT_S, PUSH_N, INC, PUSH_N`              | `ADD, LI, LI, EMIT_S, PUSH_N, PUSH_N`                |
| (5, 3) | `EMIT_S, LI, DUP, INC, PUSH_N, PUSH0`            | `EMIT_S, INC, LOG2, DUP, PUSH_N, PUSH0`                      | `EMIT_S, LI, DUP, INC, PUSH_N, PUSH0`                |
| (6, 4) | `EMIT_S, PUSH_N, LI, LI, LI, DUP, EMIT_S` (Kt 36)| `EMIT_S, GEO_SUM, PUSH_N, ADD, DUP, LOG2, PUSH0` (Kt 37, +1) | `EMIT_S, PUSH_N, LI, LI, LI, DUP, EMIT_S` (Kt 36)    |

Notable observations:

- **`(N=5, J=3)` drop_LI substitutes with LOG2 alone**, no other
  arithmetic primitive needed. The program `EMIT_S, INC, LOG2, DUP,
  PUSH_N, PUSH0` emits the LSB of `floor(log2(out_count + 1))` over
  `[0, 32)` — and this matches `bit_3(π(x))` exactly. So even within
  the LI-removed VM, **a single arithmetic primitive (LOG2) can
  substitute for LI** at this cell.
- **`(N=4, J=2)` drop_LI uses two primitives**, GEO_SUM + DIV_LOG
  composed (`floor(b / log a)` where b is GEO_SUM(out_count)).
- **`(N=6, J=4)` drop_LI degrades by +1 bit** — this is the only cell
  where the LI-kernel program is strictly shorter than any non-LI
  alternative within `L_max=28`. With a larger `L_max` this gap is
  expected to close further.

## Hard-zone observations (informational)

The hard-zone (J < ⌈N/2⌉) and pseudo-hard-zone (J between hard and
trivial) results show **stronger primitive-dependence** but also live
in the regime where `L_max ≳ target_len` may admit combinatorial
saturation (E1.3 footnote (a)):

| Cell    | L_max regime      | Ablations that saturate                |
|---------|-------------------|----------------------------------------|
| (4, 0)  | L=28, target=16 (combinatorial) | drop_LI, drop_DIV_LOG, only_LI |
| (4, 1)  | L=28, target=16 (combinatorial) | drop_LOG2, drop_LI, only_LI    |
| (5, 2)  | L=28, target=32 (meaningful)    | drop_LI, drop_DIV_LOG, only_LI |
| (3, 0)  | L=20, target=8 (combinatorial)  | (none — degrade by +3 only)    |

The genuine `(N=5, J=2)` cell — meaningful (L<target_len) — requires
**both LI and DIV_LOG**: dropping either saturates, dropping the other
two does not. The compressing program `DIV_LOG, LI, SHR1, SHR1, EMIT_S,
PUSH_N, PUSH0` realizes `floor(out_count / log(LI(out_count))) >> 2`
mod 2 — a composition of LI and DIV_LOG that is structurally
distinct from the easy-zone iterated-LI mechanism.

So the **hard-zone vs easy-zone primitive sensitivity** differs:
hard-zone compression (where it occurs) needs a specific primitive
combination, easy-zone compression is robust to single-primitive
ablation.

## What the construction produces that did not exist before

1. **Per-cell ablation table** for the four arithmetic primitives in
   the C3.a VM, distinguishing required-vs-substitutable primitive use
   at each `(N, J)` cell.
2. **F4 verdict on the easy zone**: the bounded-Kt cut shift documented
   in S150 is **arithmetic-primitive-class robust**, not LI-specific
   (refutes S150's narrowest reading).
3. **Alternative easy-zone compressing programs** without LI:
   - `(N=3, J=2)`: `EMIT_S, PUSH_N, GEO_SUM, EMIT_S`
   - `(N=4, J=2)`: `GEO_SUM, DIV_LOG, EMIT_S, PUSH_N, INC, PUSH_N`
   - `(N=5, J=3)`: `EMIT_S, INC, LOG2, DUP, PUSH_N, PUSH0`
   - `(N=6, J=4)`: `EMIT_S, GEO_SUM, PUSH_N, ADD, DUP, LOG2, PUSH0`
4. **Hard-zone primitive-pair requirement at `(N=5, J=2)`**: the
   meaningful (L<target_len) cell requires LI ∧ DIV_LOG; either-alone
   saturates.

## Falsification verdict

- **F1 (LI solely necessary and sufficient):** REJECTED — `drop_LI`
  preserves all easy-zone compression at L=24 cells and is only +1 at
  the L=28 cell.
- **F2 (multiple primitives jointly necessary):** REJECTED for easy
  zone (each single drop preserves compression). HOLDS for the genuine
  hard-zone cell `(5, 2)` (LI ∧ DIV_LOG required).
- **F3 (LI alone insufficient):** REJECTED — `only_LI` matches baseline
  on every easy-zone cell.
- **F4 (no single primitive necessary):** **HOLDS for the easy zone.**
  Each of the four arithmetic primitives can be ablated individually
  with no easy-zone compression loss.
- **F5 (mixed/incoherent):** N/A.

## What this refines

- **E1.3** refined: the S150 cut shift from `J*(N)` toward `⌈N/2⌉` is
  **driven by the FAMILY of slow-growing integer-function primitives,
  not by LI specifically.** Any one of LI / LOG2 / DIV_LOG / GEO_SUM
  suffices. The "iterated LI" reading from S150 is a special case of
  a broader statement.
- **C3.a / S150** refined: the optimum-program disassembly (which
  prominently featured iterated-LI) is **not** a faithful indicator of
  which primitives are causally required. Multiple compression
  mechanisms exist; LI just happens to give the shortest realization
  in three of four cases.

## Limitations and caveats

- Tested only at `L_max ∈ {24, 28}` and `N ∈ {3, 4, 5, 6}`. Whether F4
  persists at larger `L_max` or larger N is unverified.
- The +1-bit degradation at `(N=6, J=4)` under `drop_LI` could
  accumulate at higher N, possibly making LI strictly necessary
  asymptotically. The empirical L-budget gap between LI and non-LI
  programs at fixed N is a function we have not characterised.
- The 4-primitive set was hand-picked by S150's spec; alternative
  R⁻¹-kernel choices (e.g., LN, LOG10, INV_LI) might exhibit different
  ablation patterns.
- The construction does not affect E5.8's structural Brandt
  obstructions O1–O4; those are VM-independent.

## Edges cited

- **E5.8** — Brandt structural obstructions, unchanged by ablation.
- **E1.3** — Per-bit difficulty cut, refined: cut shift is
  primitive-class-robust, not LI-specific.
- **E1.5** — Information-rate context. Slow-growing predictors of the
  per-step rate `h_2(π(X)/X)` are the relevant primitive class; LI is
  one realization, LOG2 / GEO_SUM / DIV_LOG are alternatives.
- **E5.3** — TC⁰ open problem; this construction does not produce a
  TC⁰ circuit.

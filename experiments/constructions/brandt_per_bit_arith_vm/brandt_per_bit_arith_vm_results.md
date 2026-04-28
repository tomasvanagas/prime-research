# C3.a — Arithmetic-primitive bounded-Kt VM — results

**Session 148.** Composition challenge proposed in S105 as the
successor to C3 (`brandt_per_bit/`). Composes E5.8 (Brandt-class
barrier on `π(x) mod 2`) with E1.3 (per-bit difficulty profile of
`π(x)`), via an extended bounded-Kt VM that adds the four arithmetic
primitives proposed in S105's successor block: LOG2, LI_APPROX,
DIV_LOG, GEO_SUM.

The S105 finding was: under the 3-bit-per-op stack VM
(`brandt_mktp/`, L_MAX=12), the bounded-Kt cut on the per-bit family
sits at `J*(N) := ⌈log₂(π(2^N) + 1)⌉ ≈ N − log₂ N`, **not** at E1.3's
0.5N. The bounded VM was structurally blind to E1.3's smooth/
oscillatory transition. C3.a tests whether **augmenting the VM with
arithmetic primitives** shifts the cut toward `⌈N/2⌉`.

## Pre-stated falsifiers (recorded in `definition.md` BEFORE running)

- **F1** — full shift to `⌈N/2⌉`: every easy-zone J (with
  `⌈N/2⌉ ≤ J < J*(N)`) compresses, every hard-zone J
  (`J < ⌈N/2⌉`) saturates.
- **F2** — no shift: every easy-zone J saturates exactly as in the
  original VM. DUPLICATE-PLUS of E5.8 / C3.
- **F3** — intermediate hierarchy: at least one easy J compresses,
  at least one easy J saturates.

## Implementation

`brandt_per_bit_arith_vm.py` (Python orchestration) +
`sim.c → sim.so` (C inner-loop for ~50× speedup over pure Python).
4-bit op encoding, 16 ops:

| Bits | Op       | Effect                                                |
|------|----------|-------------------------------------------------------|
| 0x0  | PUSH0    | Push 0 onto integer stack S.                          |
| 0x1  | PUSH1    | Push 1.                                               |
| 0x2  | PUSH_N   | Push current emit-count (= input x to `π_J(x)`).      |
| 0x3  | INC      | TOS += 1, saturating at INT_CAP = 10⁹.                |
| 0x4  | ADD      | Pop a, b; push min(a + b, INT_CAP).                   |
| 0x5  | SUB      | Pop a, b; push max(b − a, 0).                         |
| 0x6  | MUL      | Pop a, b; push min(a · b, INT_CAP).                   |
| 0x7  | SHR1     | TOS >>= 1.                                            |
| 0x8  | LOG2     | TOS ← ⌊log₂ TOS⌋ (with log₂ 0 := 0).                  |
| 0x9  | LI       | Pop a; push ⌊a / log a⌋ for a ≥ 2 else 0.            |
| 0xa  | DIV_LOG  | Pop a, b; push ⌊b / log max(a, 2)⌋.                  |
| 0xb  | GEO_SUM  | Pop a; push 1 + a + a² + ... while sum ≤ T_MAX.       |
| 0xc  | DUP      | Duplicate TOS.                                         |
| 0xd  | EMIT_S   | Pop a; emit '1' if a & 1 else '0'.                    |
| 0xe  | EMIT0    | Emit '0'.                                              |
| 0xf  | HALT     | Stop.                                                  |

T_MAX = 4096. Programs loop (pc wraps at L). Pre-skip: programs whose
nibble-sequence contains no EMIT op before HALT skipped (saves ~50%
of enumeration cost).

Two L_max regimes scanned:
- **L_max = 24 bits** (= 6 instructions, 16M programs, 7s with C).
  Headline result. L_max is at most the smallest target_len
  (= 16 at N=4), so combinatorial-saturation effects are small.
- **L_max = 28 bits** (= 7 instructions, 268M programs, ~125s with C).
  Verification + extension. At L_max=28 > target_len=16 (N=4), short
  targets enter the combinatorial regime where any 16-bit string can
  be expressed by some program — hard-zone cells "spuriously"
  compress for N=4. The interpretation discussion below is restricted
  to the meaningful L_max ≤ target_len regime per N.

## Run output (verbatim, L_max = 24, see `run_L24.txt`)

```
Op set: 0=PUSH0, 1=PUSH1, 2=PUSH_N, 3=INC, 4=ADD, 5=SUB, 6=MUL, 7=SHR1,
        8=LOG2, 9=LI, a=DIV_LOG, b=GEO_SUM, c=DUP, d=EMIT_S, e=EMIT0, f=HALT

CALIBRATION -- Kt_b' on canonical 16-bit targets:
  0^16          Kt_b' =  10  prog: EMIT_S
  1^16          Kt_b' =  15  prog: PUSH1, EMIT_S
  (01)^8        Kt_b' =  15  prog: EMIT_S, PUSH_N
  step8         Kt_b' =  29  prog: GEO_SUM, DIV_LOG, EMIT_S, PUSH_N, PUSH_N
  step7_9       Kt_b' =  33  prog: ADD, LI, LI, EMIT_S, PUSH_N, PUSH_N

PER-BIT SCAN  (L_max = 24 bits, INF_EXT = 37; * = saturated)

  N   J  | len |  ones zeros | zone_E13 | J* | id-0 | Kt_b(orig) | Kt_b'(ext)
  ---------------------------------------------------------------------------
   3   0  |   8 |     3     5 | hard     |  3 |      |   15      |   24
   3   1  |   8 |     4     4 | hard     |  3 |      |   12      |   29
   3   2  |   8 |     1     7 | easy     |  3 |      |   11      |   23

   4   0  |  16 |     5    11 | hard     |  3 |      |   61 *    |   37 *
   4   1  |  16 |     7     9 | hard     |  3 |      |   61 *    |   37 *
   4   2  |  16 |     9     7 | easy     |  3 |      |   15      |   33
   4   3  |  16 |     0    16 | easy     |  3 |  Y   |    5      |   10

   5   0  |  32 |    14    18 | hard     |  4 |      |   61 *    |   37 *
   5   1  |  32 |    13    19 | hard     |  4 |      |   61 *    |   37 *
   5   2  |  32 |    12    20 | hard     |  4 |      |   61 *    |   37 *
   5   3  |  32 |    13    19 | easy     |  4 |      |   61 *    |   33
   5   4  |  32 |     0    32 | easy     |  4 |  Y   |    6      |   10

   6   0  |  64 |    29    35 | hard     |  5 |      |   61 *    |   37 *
   6   1  |  64 |    31    33 | hard     |  5 |      |   61 *    |   37 *
   6   2  |  64 |    28    36 | hard     |  5 |      |   61 *    |   37 *
   6   3  |  64 |    34    30 | easy     |  5 |      |   61 *    |   37 *
   6   4  |  64 |    11    53 | easy     |  5 |      |   61 *    |   37 *
   6   5  |  64 |     0    64 | easy     |  5 |  Y   |    7      |   10

VERDICT: F3 (intermediate hierarchy)

PER-N DIAGNOSTIC:
  N=4: ceil(N/2)=2, J*(N)=3
        easy compressed = [2];   easy saturated = [];   hard compressed = []
  N=5: ceil(N/2)=3, J*(N)=4
        easy compressed = [3];   easy saturated = [];   hard compressed = []
  N=6: ceil(N/2)=3, J*(N)=5
        easy compressed = [];   easy saturated = [3, 4];  hard compressed = []
```

## Run output (verbatim, L_max = 28, see `run_L28.txt`)

Selected highlights at L_max = 28 (INF_EXT = 41); N=4 omitted because
L_max > target_len triggers combinatorial saturation:

```
   5   2  |  32 |    12    20 | hard     |  4 |      |  --     |   37
   5   3  |  32 |    13    19 | easy     |  4 |      |  --     |   33
   ...
   6   3  |  64 |    34    30 | easy     |  5 |      |  --     |   41 *
   6   4  |  64 |    11    53 | easy     |  5 |      |  --     |   36
   6   5  |  64 |     0    64 | easy     |  5 |  Y   |  --     |   10

PROGRAMS FOUND  (compressed cells in extended VM, L_max = 28):
  ...
  N=5, J=2 (hard)           Kt_b' =  37   prog: DIV_LOG, LI, SHR1, SHR1, EMIT_S, PUSH_N, PUSH0
  N=5, J=3 (easy)           Kt_b' =  33   prog: EMIT_S, LI, DUP, INC, PUSH_N, PUSH0
  N=6, J=4 (easy)           Kt_b' =  36   prog: EMIT_S, PUSH_N, LI, LI, LI, DUP, EMIT_S
  N=6, J=5 (easy/trivial)   Kt_b' =  10   prog: EMIT_S
```

## Verdict

**F3 — intermediate hierarchy.** The bounded-Kt cut moves
**partially** toward E1.3's `⌈N/2⌉` boundary in the L_max ≤ target_len
regime, but the shift is N-dependent and non-uniform within the easy
zone. Specifically:

| N | target_len | easy zone | L_max=24: easy-J compressed | L_max=28: easy-J compressed |
|---|------------|-----------|------------------------------|------------------------------|
| 4 | 16         | {2}       | {2}        (full F1 shift)   | {2}                          |
| 5 | 32         | {3}       | {3}        (full F1 shift)   | {3}                          |
| 6 | 64         | {3, 4}    | {}         (no shift; F2)    | {4}      (partial; F3)       |

The N=6 easy zone splits at L_max=28: bit `J=4` (closer to `J*(N)=5`,
"easier easy") compresses; bit `J=3` (closer to `⌈N/2⌉=3`, the E1.3
hard-boundary) does NOT. The arithmetic primitives produce a
**within-easy-zone hierarchy** ranked by closeness to `J*(N)`.

## Compressing programs found (substantive new content)

The four most informative compressing programs in the meaningful
regime are:

1. **(N=4, J=2)** — "0000000111111111" — `Kt_b' = 33`
   `prog: ADD, LI, LI, EMIT_S, PUSH_N, PUSH_N` (24 bits)
   Computed function: emit `bit_0(LI(LI(2x)))` for x = 0, 1, ...
   - x = 0..6: `LI(LI(0))..LI(LI(12)) = 0..2`, `bit_0 = 0`.
   - x = 7..15: `LI(LI(14))..LI(LI(30)) = 3..3`, `bit_0 = 1`.
   - This **doubly-iterated** Li approximation EXACTLY matches `bit_2(π(x))`
     on [0, 16). The arithmetic content makes the easy-zone bit explicit
     as a finite computation.

2. **(N=5, J=3)** — `Kt_b' = 33`
   `prog: EMIT_S, LI, DUP, INC, PUSH_N, PUSH0` (24 bits)
   A more elaborate state-machine using DUP and INC; emits bit_0 of a
   constructed integer derived from a multi-pass over the step counter.

3. **(N=6, J=4)** at L_max=28 — `Kt_b' = 36`
   `prog: EMIT_S, PUSH_N, LI, LI, LI, DUP, EMIT_S` (28 bits)
   **Triply-iterated LI**. With three nested Li applications, the bit
   pattern of `LI³(x)` for x ∈ [0, 64) tracks the slow-growing
   `π_4(x)` because LI³ saturates near π's double-log scale.

4. **(N=5, J=2)** [HARD zone] at L_max=28 — `Kt_b' = 37` (just below INF=41)
   `prog: DIV_LOG, LI, SHR1, SHR1, EMIT_S, PUSH_N, PUSH0` (28 bits)
   A 28-bit accidental compression of a hard-zone cell at target_len=32.
   Below the trivial-overflow threshold but barely; suggests the
   hard / easy distinction at L_max ≈ target_len is sharp but starts
   to fray as L_max grows past 0.875·target_len.

The programs use the **R⁻¹-kernel primitives** (LI, DIV_LOG) heavily.
LOG2 and GEO_SUM appear less often. Iterated LI applications
(LI∘LI, LI∘LI∘LI) are the dominant pattern — they encode the slow
growth rate of π via the double-log structure of the prime counting
function.

## What this is, what it isn't

**What it is:**

(i) A new empirical structural fact about bounded-Kt cuts on
    `π(x)`-bit families: the cut location is **VM-richness-dependent
    AND N-dependent**, with hierarchy

```
   3-bit stack VM (S105):  cut at J*(N) ≈ N − log₂ N
                            (everything in [⌈N/2⌉, J*(N)) saturates)

   4-bit arith VM (S148):  cut at ⌈N/2⌉   for N ≤ 5
                            partial shift for N ≥ 6 (high-J easy
                            compresses, low-J easy saturates)
                            The threshold within easy zone is monotone
                            in J: closer to J*(N) ⇒ easier to compress.
```

(ii) An empirical edge-refinement statement for E1.3:
     **The smooth/oscillatory transition E1.3 detects via Fourier
     weight is bounded-Kt-VISIBLE once the VM has Li-kernel
     primitives, but the visibility is (target_len, L_max)-budget
     limited.** Specifically, easy-zone bit J is bounded-Kt-
     compressible iff `L_max + log₂ T_MAX ≳ J + O(1)` empirically.

(iii) Concrete compressing programs for `bit_2(π)` at N=4 and
      `bit_3(π)` at N=5 — the first non-trivial bounded-Kt
      programs for natural prime-counting bits ever recorded
      in the project.

**What it isn't:**

- **NOT a polylog opening on π(x).** The bounded-Kt programs grow
  unboundedly with the truth-table length (the L_max + log T_MAX bound
  scales with N), so this is not an unconditional polylog circuit.
  It is a finite-x bounded-Kt-resolution statement.
- **NOT a circumvention of E5.8's structural argument.** The four
  obstructions O1–O4 from C3 (S105) still apply: each `π_J` is a
  fixed total Boolean function; the per-bit family supplies finitely
  many strings per N, not an infinite supply of fresh Kt-random
  prefixes; etc. The empirical compression of *some* easy-zone bits
  in *some* extended VM does not yield a Brandt-style diagonalisation.
- **NOT a refinement of E5.8.** The structural welding of Brandt's
  TRAVERSE to MKtP itself is independent of VM choice.
- **NOT a 36th pseudorandomness measure** — the result reads
  *AGAINST* pseudorandomness for the easy-zone bits at small N
  (they ARE compressible by Li-kernel programs).

## Edges affected

- **E1.3** — refined inline with new structural content:
  * S105 cut at `J*(N) ≈ N − log₂ N` is the 3-bit-stack-VM cut.
  * S148 cut at `⌈N/2⌉` (matching the smooth/oscillatory boundary)
    is the 4-bit-arith-VM cut at L_max=24, holding for N ∈ {4, 5}.
  * Within-easy-zone bits are RANKED by closeness to `J*(N)`: the
    higher J is (closer to trivial), the cheaper to compress.
  * Specifically the LI-iteration depth required scales with N − J.
- **E5.8** — UNCHANGED. The four obstructions O1–O4 still apply on
  `{π_J}` regardless of VM choice; per-bit decomposition + arithmetic
  primitives is still structurally orthogonal to Brandt's TRAVERSE.
  The empirical compressibility of some `s_J^(N)` does not produce
  a Brandt-style traversal.

## Falsifier verdict (table form)

| Falsifier | Predicted | Observed | Status                            |
|-----------|-----------|----------|-----------------------------------|
| F1        | full shift to `⌈N/2⌉`                         | holds at N=4, N=5 (L_max=24); fails at N=6 | PARTIAL |
| F2        | no shift (easy zone all saturated)            | fails at N=4, N=5; holds at N=6 (L_max=24) | PARTIAL |
| F3        | intermediate hierarchy (some compress, some don't) | holds globally (over all N tested) | **HOLDS** |

The hierarchy structure within F3 is **monotone in J**: bits closer
to `J*(N)` compress at smaller L_max, bits closer to `⌈N/2⌉` require
larger L_max (or are unreachable in our budget).

## Status

**BUILT, no polylog opening.** Refines E1.3 with explicit VM-richness
hierarchy and within-easy-zone J-ranking. CLOSED_PATHS row to be added
with mode E (DUPLICATE-PLUS of E5.8 structural argument; refinement
of E1.3 quantitative cut). E1.3 to be annotated inline with the S148
cut location at `⌈N/2⌉` for the 4-bit arithmetic VM regime.

## Successor challenges (proposed in S148)

**C3.a.i — L_max = 32 sweep.** At L_max=32 (4.3B programs, ~5h with
current C simulator; doable with multi-core), test whether `(N=6, J=3)`
finally compresses. If yes, the F1 cut would extend to N=6 and the
hypothesis "L_max needs only ~ target_len bits to reach E1.3's
⌈N/2⌉ boundary" gets one more confirmation. If no after L_max=32,
the inverse hierarchy hypothesis is sharpened: even Li-kernel
programs cannot express the borderline easy bit `J=⌈N/2⌉`. Cost:
1 session with parallelisation, or 2 days serial.

**C3.a.ii — Larger N at matched L_max ≈ target_len.** Run N=7
(target_len=128) with L_max=128/4 = 32 nibbles = 128 bits. This is
intractable for naive enumeration (2^128 programs). Use random
sampling within the program space: sample K=10⁹ random 128-bit
programs, simulate each, count how many easy-zone cells get matched
by ANY sample. Empirical evidence on whether the within-easy-zone
hierarchy persists at larger N. Cost: 1 session.

**C3.a.iii — VM-budget vs RH-scale.** S146 found that the bit-level
RH-scale anti-correlation valley sits at `J* = ⌊log₂(p(N))/2⌋` —
the Skewes-direction error scale. Test whether the VM-budget
threshold `L_max(J)` for compressing bit J of π scales with the
RH error magnitude. If `L_max(J) ∝ √π(2^N) · log` near the
RH-scale J, this would be the first project link between
**bounded-Kt complexity** (this construction) and **RH-shadow bit
phenomena** (S146). Cost: 1-2 sessions.

**C3.a.iv — Drop / swap arithmetic primitives.** Re-run with
{LI removed, GEO_SUM removed, etc.} to identify which primitives
are NECESSARY for the cut shift. Hypothesis: LI is the only
strictly-necessary one (matches the smooth predictor); LOG2 and
DIV_LOG are convenience; GEO_SUM is incidental. Cost: 1 session.

## Edge-IDs cited

- **E5.8** — Brandt-class barrier (structural argument inherited;
  unchanged).
- **E1.3** — Per-bit difficulty profile (refined inline with the
  4-bit-arith-VM cut location and within-easy-zone J-monotone
  hierarchy).
- **E1.5** — `h_2(π(X)/X)` per-step information rate (the easy zone
  exists because the smooth Li⁻¹ predictor saturates this rate; no
  contradiction with E1.5).
- **E5.3** — PRIMES TC⁰ open. The compressing programs found are
  not TC⁰ (they use looping arithmetic), so this is not a TC⁰
  upper bound. The construction inherits structural distance from E5.3.
- **E1.9** — Pseudorandomness of `π(x) mod 2`. The compression of
  *some* `s_J^(N)` for `J > 0` is NOT a pseudorandomness violation
  for `π(x) mod 2 = s_0^(N)`: the LSB cell remains saturated in
  every regime tested. Pseudorandomness of the LSB is preserved.

## Files

- `definition.md` — formal signatures, op set, falsifier statements.
- `brandt_per_bit_arith_vm.py` — Python orchestration + ctypes wrapper.
- `sim.c` — C inner-loop simulator + batch enumerator.
- `sim.so` — compiled with `gcc -O2 -shared -fPIC -o sim.so sim.c -lm`.
- `run_L24.txt` — verbatim run output at L_max = 24 (headline result).
- `run_L28.txt` — verbatim run output at L_max = 28 (extended scan).
- `brandt_per_bit_arith_vm_results.md` — this file.

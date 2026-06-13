# pcs_commit.py — tensor multilinear PCS with a sub-√x verifier (S508, S514)

> **S514 update — field-size-free row code.** The RS row code needs the codeword
> length `N = blowup·k ≤ q`, a demo constraint that capped the committable table
> size and forced the large runs onto the 61-bit BIG_Q. `pcs_commit.py` now offers
> a **Brakedown-style linear-time expander code** behind the SAME
> `commit/prove/verify` interface, gated by `code=` (`"rs"` default keeps the RS
> path verbatim; `"expander"` removes the `N≤q` requirement). Same sub-√x verifier
> slope, same 4-class cheat panel, over q & BIG_Q & SMALL_Q, and a concrete tiny-
> field demo (q=17) where RS *cannot* commit but the expander does. See the S514
> section at the bottom.

## What this is

A hash-based multilinear polynomial commitment (tensor / Ligero–Brakedown
style) whose **evaluation-proof verifier is sub-√x** — `O(√(2^nb)·polylog) =
O(x^{1/4}·polylog)`. It exists to discharge the chain's **last** verifier work
whose cost grows with `x`: the two *one-time* `S_0` base closes in
`compressed_layer.run_chain`, each a direct `mle_eval` over a `2^nb = O(√x)`-size
table (`C == mle_eval(S0, z)`, `compressed_layer.py:1047`).

After S506 (`batch_ub`) the per-layer verifier is leaf-eval-free; S507 measured
the whole-chain verifier op-count at `Θ(√x)` (per-layer 0, one-time `2·2^nb`).
This module collapses that one-time `2^nb` term to `O(x^{1/4})`, so the
**whole-chain verifier is sub-√x (succinct under a hash assumption)**.

## Construction (tensor / Ligero–Brakedown)

For a multilinear `S~` over `nb` vars, table `S` of `2^nb` values (MSB-first):

- Reshape `S` into an `r×k` matrix `M` (`r=2^{n1}`, `k=2^{n2}`, `n1+n2=nb`, high
  bits = rows, low bits = cols → `M = S.reshape(r,k)`).
- `eq~(pt,w)` factors over the split, so with `a = eq_table(pt_hi)` (len `r`),
  `b = eq_table(pt_lo)` (len `k`):  `S~(pt) = a^T M b`.
- **commit**: Reed–Solomon-encode each *row* of `M` (degree-`<k` message → `N =
  blowup·k` evaluation points), giving `Mhat (r×N)`; Merkle-commit the `N`
  *columns* (sha256 = the CRH/random-oracle stand-in). Commitment = the root.
- **prove** (Fiat–Shamir): `rho ← RO(root,pt,claimed)` (len `r`); send the two
  combined messages `v = a^T M` (evaluation) and `w = rho^T M` (proximity), then
  open `t` columns `Q ← RO(root,pt,claimed,v,w)` with Merkle paths.
- **verify**: (1) evaluation binding `⟨v,b⟩ == claimed`; (2) per queried column
  `c`: Merkle membership of `Mhat[:,c]` to `root`, AND `Enc(v)[c] == ⟨a,Mhat[:,c]⟩`
  AND `Enc(w)[c] == ⟨rho,Mhat[:,c]⟩` (encoding is linear, so for an honest
  codeword matrix `Enc(a^T M)[c] = a^T Mhat[:,c]`).

Verifier touches `O(k)` (read `v,w`, `⟨v,b⟩`) + per column two len-`r` dot
products + two len-`k` encodings = `O(t·(r+k))` field ops. With `r,k ~ 2^{nb/2} ~
x^{1/4}` and `t` fixed → `O(x^{1/4}·polylog)`. The full `2^nb` table never touches
the verifier (it lives only in the prover's commit/open).

## Measured — `python3 pcs_commit.py --bench` (default q = 2³¹−1)

```
 nb   2^nb     r     k    N  commit_vops  direct_2^nb  ratio  proof_cols
  4     16     4     4   16          260           16  0.06x        16
  6     64     8     8   32         1032           64  0.06x        32
  8    256    16    16   64         2064          256  0.12x        32
 10   1024    32    32  128         4128         1024  0.25x        32
 12   4096    64    64  256         8256         4096  0.50x        32
 14  16384   128   128  512        16512        16384  0.99x        32

fitted slope d log2(ops)/d nb (last 4 pts):  commit = 0.500   direct = 1.000
in x: alpha_commit ~ 0.250 (target ~0.25),   alpha_direct ~ 0.500 (~0.5).
```

**The headline is the SLOPE.** `commit_vops` grows as `Θ(2^{nb/2}) = Θ(x^{1/4})`
(slope 0.5 per nb = 0.25 per `n`), the direct close as `Θ(2^nb) = Θ(√x)`. Same
slopes over `--field big` (BIG_Q = 2⁶¹−1).

**Honest crossover.** The *absolute* op count crosses the direct close only at
`nb ≈ 14` (`commit_vops 16512` vs `direct 16384`), i.e. `n ≈ 28`. Below that the
commit verifier does MORE work (the `t`-query constant dominates): the win is
asymptotic, not at reachable demo `n` — exactly the S506 wall-clock pattern. The
falsifiable claim is the slope (sub-√x), which holds at every measured `nb`.

## Soundness — cheating-prover tests (selftest, over q, BIG_Q, SMALL_Q; nb ≤ 8)

All rejected:
- **wrong claimed value** — fails the evaluation binding `⟨v,b⟩ == claimed`
  (the honest `v` sums to the true value).
- **forged opening** — a `v'` bent at one coordinate so `⟨v',b⟩ = wrong` (passes
  binding by construction, asserted) is caught by the column-consistency checks
  (`Enc(v')[c] ≠ ⟨a,Mhat[:,c]⟩`; two degree-`<k` polys agree on `<k` of `N`
  columns, reject prob `≥ 1−((k−1)/N)^t`).
- **tampered (non-codeword) committed row** — Merkle root rebuilt to a valid
  commitment of garbage; caught by the random proximity columns.
- **tampered revealed column** (single entry) — fails its Merkle path.

Honest openings agree with `mle_eval` at every tested `nb`; the commitment is the
same map as the chain's `pad`+`mle_eval` base close (selftest case 7).

## Integration — `compressed_layer.run_chain(commit_base=True)`

The base close branches: default = direct `mle_eval` (unconditional, `O(√x)`);
`commit_base` = two tensor-PCS opens (verifier sub-√x, tallied in
`stats['vcommit_ops']`; proof comm added to `comm`). Verdict UNCHANGED — honest
`claimed π == sieve`, `delta_pi` and self-consistent liar rejected — over q AND
BIG_Q (selftest case 23, `n ≤ 12`).

`python3 compressed_layer.py --bench-verifier-ops` (config **(d)** added):

```
config         alpha_ops   expectation
a no-pcs          0.961     ops~x      (per-layer leaf closes)
b pcs             0.998     ops~x      (Ub openings survive)
c pcs+bt+ub       0.500     ops~sqrt x (per-layer leaf-eval-free; one-time 2·2^nb)
d +commit         0.258     ops~x^1/4  (one-time base via tensor PCS)  ← SUB-√x
```

The whole-chain verifier leading term drops `Θ(x) → Θ(√x) → Θ(x^{1/4})`. `comm`
slope stays `~0.60` (Õ(√x)) — the commit proof adds only an `O(x^{1/4})` term.

## Honest scope / caveats

- **Computational, not unconditional.** This trades the otherwise-unconditional
  Õ(√x) verifier for full succinctness under a collision-resistant-hash /
  random-oracle assumption (sha256 stand-in). The rest of the chain verifier is
  unconditional.
- **Soundness is the standard Ligero/Brakedown argument**, verified here
  *empirically* by the cheating-prover tests — not a machine-checked proof.
- **RS field-size — now removed by `code="expander"` (S514).** RS needs `N ≤ q`
  distinct points (fine at demo `n`, but caps the table size and forces large runs
  onto a 61-bit field). The field-size-free Brakedown-style expander code (S514
  section below) plugs in behind the SAME `commit/open/verify` interface; the
  tensor reduction and the `O(x^{1/4})` verifier cost are unchanged.
- **Leading-term metric.** `vops`/`vcommit_ops` counts the verifier's dominant
  field-op loops (`O(t(r+k))`); the `O(t·log N)` Merkle hashing is also sub-√x.

## What would falsify this

`verify` disagreeing with `mle_eval` on an honest opening; a forged opening /
wrong claimed value / tampered codeword row accepted above the field bound; the
verifier op-count slope `→ 0.5` instead of `→ 0.25` (i.e. not sub-√x); or
`run_chain(commit_base=True)` changing the chain verdict / `claimed π`.

---

## S514 — the field-size-free Brakedown-style expander row code (`code="expander"`)

### Why

The RS row code requires the codeword length `N = blowup·k ≤ q` (every codeword
symbol is a distinct field point). At demo `nb` this is fine, but it **caps the
committable table size** and, in `large_x_benchmark.py`, is one reason the runs
use the 61-bit `BIG_Q`. A production-scale commitment (`n ≳ 120`, `nb ≳ 60`) over
a small or fixed field needs a code whose length is independent of `q`.

### What

A **systematic, linear, linear-time-encodable** expander code (Spielman / Druk–
Ishai recursive construction, as instantiated in Brakedown, Golovnev–Lee–Setty–
Thaler–Wahby 2021), plugged in behind the **same** `commit/prove/verify`:

```
Enc_n(x) = x  ||  Enc_{nm}(A·x)  ||  B·Enc_{nm}(A·x),   nm = ceil(n/2)
```

with `A (nm×n)`, `B (pn×ne)` **column-regular** sparse matrices over `F_q`
(`coldeg=4` nonzeros per column) and a dense systematic base case `x || R·x` for
`n ≤ 4`. The matrices are derived deterministically from a fixed public seed +
the size `n` (column indices, `q`-free) with values reduced mod `q`, so prover
(commit) and verifier (re-encode) build the **identical public linear map**.

The only code-specific code paths are **commit** (encode each row) and **verify**
(re-encode the combined messages `v,w` once, `O(N)=O(k)`, then index the queried
columns). The tensor identity `S~(pt)=a^T M b` and the homomorphic column check
`⟨a,Mhat[:,c]⟩ = Enc(a^T M)[c] = Enc(v)[c]` hold for **any** linear `Enc`, so the
soundness argument and the `O(t(r+k))=O(x^{1/4})` verifier are unchanged.

**Column-regularity is load-bearing.** A plain low-density generator matrix (each
output row picks random input columns) leaves some input coordinate unused → a
basis codeword of weight 1 → the forge cheat at that coordinate is essentially
uncatchable. Making every *input column* hit `coldeg` output rows guarantees every
input influences the parity. Measured min basis-codeword relative weight (the
quantity that governs the forge cheat — the forged-`v` difference is
`δ·Enc(e_{j0})`):

```
              row-regular (cols missed)   column-regular (used here)
  k=128       rel ~ 0.003                 rel ~ 0.45–0.47
```

### Measured — `python3 pcs_commit.py --bench --code expander`

```
 nb   2^nb     r     k    N  commit_vops  direct_2^nb  ratio  proof_cols
  4     16     4     4    8          100           16  0.16x         8
  6     64     8     8   20          488           64  0.13x        20
  8    256    16    16   44         1488          256  0.17x        32
 10   1024    32    32   92         3136         1024  0.33x        32
 12   4096    64    64  188         6464         4096  0.63x        32
 14  16384   128   128  380        13152        16384  1.25x        32

fitted slope d log2(ops)/d nb (last 4 pts):  commit = 0.524   direct = 1.000
in x: alpha_commit ~ 0.262 (target ~0.25),   alpha_direct ~ 0.500 (~0.5).
```

- **Slope sub-√x (0.524, α≈0.262)** — same as RS (0.500); the expander keeps the
  `Θ(x^{1/4})` verifier. Identical slope over q, BIG_Q, and SMALL_Q (the op count
  is field-independent — only the indices, not the values, set the structure).
- The expander codeword is **shorter** than RS's (e.g. `N=380` vs RS `512` at
  nb=14, blowup 4), so `commit_vops` is actually *lower* at large nb — the verifier
  re-encodes `v,w` once (`O(N)`) instead of per-column (`O(t·k)`).

### Soundness — the 4 cheat classes, over q & BIG_Q & SMALL_Q (selftest)

Same panel as RS, all rejected with `code="expander"`: **wrong claimed value**
(eval binding), **forged opening** (column consistency — the difference is
`δ·Enc(e_{j0})`, weight `≥ δN` with measured `δ≈0.45`, caught w.p. `≥1−(1−δ)^t ≈
1−7e-9` at `t=32`), **tampered codeword row** (proximity columns), **tampered
revealed column** (Merkle path). Honest openings agree with `mle_eval`; same map
as the chain base close (selftest case 7, both codes); measured min basis-codeword
relative weight `≥ 0.3` asserted for `k∈{8,16,32,64}` over q & BIG_Q
(`_exp_min_basis_rel_weight`).

### The field-size-free win — concrete (selftest case 9)

Over a **tiny prime q=17**, at `nb=8` (`k=16` → RS `Ncode = max(64,17) = 64 > 17`):

```
RS  over q=17, nb=8:  refuses to commit  (RS needs Ncode <= q distinct points)
EXP over q=17, nb=8:  commits (N=44) & verifies = True; wrong-claim & forged-
                      opening both rejected
```

RS **cannot** commit (codeword longer than the field); the expander commits,
verifies, and rejects the cheats. This is exactly the `N≤q` constraint the
production face needs removed.

### Chain integration — `run_chain(commit_base=True, commit_code="expander")`

`commit_code` selects the row code for the two `S_0` base opens (default `"rs"`
keeps every prior artifact verbatim). End-to-end, `--n 12 --delegate --structured
--field big --pcs --commit-base --commit-code expander`: honest `claimed
π(4095)=564 == sieve`, `delta_pi` + self-consistent liar + all 13 layer cheats
rejected. Selftest **case 24** asserts verdict UNCHANGED vs the RS base over q AND
BIG_Q at `n∈{8,10,12}` (per-layer still leaf-free, base opens still sub-√x:
`vcommit_ops == vleaf_ops_ot`), cheats rejected.

### Honest scope (S514)

- **Distance is the security parameter, and it is MEASURED here, not proven from
  the asymptotic expander bound.** Demo-scale parameters (`coldeg=4`, recursion
  `nm=ceil(n/2)`, parity `ceil(n/2)`) give a measured relative distance `δ≈0.45`;
  Brakedown's *analyzed* parameters give a proven `δ≈0.07` with `t≈148` queries
  for 100-bit security. The **construction shape** here is theirs; the constants
  are demo-scale and distance-measured. The soundness error of the consistency /
  proximity tests is `≤ (1−δ)^t` per the standard Ligero/Brakedown argument.
- Still computational (CRH/RO, sha256 stand-in), like the RS PCS — this swap does
  not change the unconditional-vs-computational boundary, only the field
  dependence.

### What would falsify the S514 work

The expander path disagreeing with `mle_eval` on an honest opening; any of the 4
cheats accepted with `code="expander"`; the measured min basis-codeword relative
weight collapsing toward 0 (re-introducing catchable-only-at-one-column forges);
the verifier op-count slope `→ 0.5` for the expander (not sub-√x); or
`run_chain(commit_code="expander")` changing the chain verdict / `claimed π`.

---

## S515 — distance/soundness as a measured CURVE (`--distance-sweep`, `--forge-rate`)

### Why

S514 rested the expander code's soundness on its relative minimum distance `δ`
(the forged-`v` difference is exactly `δ·Enc(e_{j0})`, so a forge false-accepts
with prob `≤ (1−δ)^t` over `t` queried columns). S514 **measured** `δ≈0.45` at a
single demo point (`k∈{8..64}`) but did not (a) confirm `δ` stays bounded below as
`k` grows toward `~512`, nor (b) check the empirical false-accept rate actually
tracks `(1−δ)^t`. This cycle turns both into measured curves, behind the
`--selftest`-gated harness.

### (a) Distance sweep — `δ = min_j |Enc(e_j)|/N` vs k (`--distance-sweep`)

`_min_basis_rel_weight(code, k, q)` (a code-agnostic generalization of S514's
`_exp_min_basis_rel_weight`) measures the min basis-codeword relative weight — the
conservative soundness parameter — for **both** row codes over q & BIG_Q & SMALL_Q,
`k` up to 512. RS basis weights via the incremental power table `c^j = c^{j−1}·c`
(O(Nk), measured); expander via `_exp_encode` of each `e_j`.

```
--- code = rs ---             --- code = expander ---
  k    N    δ (all fields)      k    N    δ (all fields)
  8   32   0.969                8   20   0.650
 16   64   0.984               16   44   0.568
 32  128   0.992               32   92   0.522
 64  256   0.996               64  188   0.473
128  512   0.998              128  380   0.447
256 1024   0.999              256  764   0.423
512 2048   1.000              512 1532   0.404
```

- **RS** `δ → 1` (the systematic monomial basis: `Enc(e_j)` vanishes only at the
  point `c=0` for `j≥1`, so `δ = (N−1)/N`). Increasing — trivially non-decaying.
- **Expander** `δ` **declines 0.650 → 0.404** over `k=8→512`, BUT the per-doubling
  decrements **shrink**: `+0.082, +0.046, +0.048, +0.026, +0.025, +0.019` (mean
  ratio `≈0.77`). A geometric-tail extrapolation puts the asymptotic floor at
  **`≈0.34`**. **δ does NOT decay to 0** — the decline decelerates toward a positive
  floor, consistent with (not a proof of) Brakedown's proven constant relative
  distance. Field-**independent** (identical over q, BIG_Q, SMALL_Q — the support
  is set by q-free indices).
- **Practical consequence** (`(1−δ)^t ≤ 2^{−100} ⇒ t ≥ 100/−log₂(1−δ)`): the query
  count for 100-bit forge soundness grows only **modestly**, `t: 67 → 134` as
  `k: 8 → 512` (RS: `t: 20 → 10`). Sub-linear in `k` precisely because `δ` floors
  above 0; a `δ→0` code would need `t→∞`.

### (b) Forge-rate Monte-Carlo — empirical accept vs `(1−δ)^t` (`--forge-rate`)

The **worst-case** adversary bends the honest `v` at the min-weight basis
coordinate `j0 = argmin_j weight(Enc(e_j))` (subject to `b[j0]≠0`), so the
difference codeword has support `δ·N`. Across many random openings (independent
Fiat–Shamir challenges via random `pt` + wrong value) we measure the fraction
accepted at each `t` — accept iff none of the first `t` FS-derived columns lands in
the support. Example (`expander`, q, nb=10, k=32, N=92, δ=0.522, 20000 trials):

```
  t   emp_accept   (1-δ)^t   hypergeom   emp/pred
  1     0.47655    0.47826    0.47826      1.00x
  2     0.22305    0.22873    0.22599      0.98x
  3     0.10195    0.10939    0.10546      0.93x
  4     0.04750    0.05232    0.04858      0.91x
  ...
 10     0.00040    0.00063    0.00034      0.64x
 12     0.00005    0.00014    0.00006      0.35x
  fitted delta_eff = 0.5616   measured delta = 0.5217
```

- The empirical accept rate **tracks `(1−δ)^t`** and stays **at or below** it
  (`emp/pred ≤ ~1.0` at every `t`), in fact hugging the tighter sampling-without-
  replacement **hypergeometric** prediction. The fitted decay rate
  `δ_eff = 0.56 ≥ δ = 0.52` — the forge dies at least as fast as the bound predicts.
- **No query correlation.** `emp/pred ≫ 1` would mean the FS columns are correlated
  (the falsifier); it is never observed (max `≈1.02`, a `~3σ` fluctuation at `t=1`).
  Holds over q AND BIG_Q (field-independent law), and for RS (`δ=0.97`, accept
  `≈0.027` at `t=1`, `0` thereafter — caught almost surely by one column).

### Soundness selftest cases (10, 11) — falsifiers

- **Case 10 (distance):** for both codes over q & BIG_Q & SMALL_Q, `k≤128`:
  `min δ ≥ floor` (RS 0.5, expander 0.40) — **δ does not decay to 0**; the
  per-doubling decline is **decelerating** (`decr_last ≤ decr_first`); expander δ is
  field-independent; `_hyper_miss` is a valid monotone `≤ (1−δ)^t` bound.
- **Case 11 (forge-rate):** worst-case Monte-Carlo (nb=6, 4000 trials, both codes):
  empirical accept `≤ (1−δ)^t` within binomial tolerance at every `t` (the
  correlation falsifier), monotone non-increasing in `t`, `t=1` accept `≈ 1−δ`, and
  fitted `δ_eff ≥ δ − 0.10` (decay not slower than predicted).

### Honest scope (S515)

- The **asymptotic constant distance is Brakedown's, proven, and cited** — it is NOT
  re-proved here. What is added is a **measurement**: at demo scale `δ` declines with
  `k` (the construction's constants, `coldeg=4`, are demo-scale, not Brakedown's
  100-bit-analyzed `δ≈0.07, t≈148`), but the decline **decelerates** and `δ` stays
  `≥0.40` to `k=512`, with a geometric-tail floor estimate `≈0.34`. The geometric
  extrapolation is an estimate, not a bound.
- The practical recipe is therefore: **measure `δ` at the actual `k`** used and pick
  `t = ⌈λ/−log₂(1−δ)⌉`. The forge-rate experiment validates that this `t` delivers
  the claimed `(1−δ)^t` soundness empirically.

### What would falsify the S515 work

`δ → 0` as `k` grows (the distance sweep showing acceleration, not deceleration, of
the decline); the empirical forge accept rate `≫ (1−δ)^t` at some `t` (FS columns
correlated); `δ_eff ≪ δ` (forge dies slower than the bound); or the law differing
across fields (it must not — the code is field-size-free).

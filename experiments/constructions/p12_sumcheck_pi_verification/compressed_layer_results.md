# compressed_layer.py — results (S493 step 2; S494 step 3)

**Step 2 of the compressed prover** (PROGRAM.md open item 1): the step-1
trace primitive (`compressed_prover_mult_trace.py`) integrated into **one real
Lucy sieve layer over the compressed O(√x)-size state**, with the
small-side/large-side dispatch the production prover needs.

**Step 3 (S494): the small-side layer + the full two-sided CHAIN.** The
large-side layer of step 2 is now joined by a small-side layer, and the two
sides are chained over all `K=π(√x)` primes to verify **π(x) end-to-end over
the compressed state** — `claimed π(x) == sieve π(x)` at every tested `x`,
with the self-consistent-liar caught at the corruption layer. **See the
"Step 3" section below; the step-2 detail that follows it is unchanged.**

**Step 4 (S495): wiring delegation — the verifier is now cleanly Õ(√x).** Both
remaining `O(nb·p)` automaton evaluations are replaced by an `O(nb·log p)`
inner protocol (`--delegate`). **See the "Step 4" section immediately below.**

**Step 5 (S496→S497): the chain wiring PROVER is now Õ(x), not Õ(x^{3/2}).** The
delegated chain prover materialized dense `2^{2l}`-size tables per chain layer
(`O(nb·p²)` per Lucy layer → `Σ_{p≤√x} nb·p² ≈ Õ(x^{3/2})`); S496 built the
structured `O(nb·p)` replacement in `lucy_dp_delegated_wiring.py`. **This cycle
(S497) propagates it into the real compressed-√x-state π(x) chain** by threading
a `structured` flag through `compressed_layer.py`'s reduction stack down to its
`inner_verify_div` / `verify_waff_value` calls (`--delegate --structured`). The
transcript is **bit-identical**, so `--selftest` asserts identical accept/reject,
comm AND claimed π(x) with structured on vs off, and the new `--prover-bench`
confirms the wiring prover flattens `p²→p`. **See the "Step 5" section below.**

**Step 6 (S499): the field lift threaded through the whole compressed chain.**
The demo prime `Q = 2³¹−1` (uint64) is sound only while every certified integer
stays below it — values up to `x` and the trace product `U·a ≈ X` — i.e. `x < Q`,
roughly `n ≲ 30`. Above that a wrap-around alias `u''·a + r'' = X + k·Q` aliases
`X mod Q` yet `u'' ≠ ⌊X/a⌋` (the S498 hole, demonstrated at the primitive level).
This cycle threads a prime-modulus parameter `q` through **every** field op of
`run_chain` and the delegated wiring (`lucy_dp_delegated_wiring.py`) — `eq_table`,
`sumcheck`, `mle_eval`, the comparator/division/affine automata, the carry-chain,
the structured prover — defaulting `q = Q` so all prior artifacts stay
**bit-identical**, with a uint64 fast path for `q ≤ 2³¹−1` and exact Python-int
**object arrays** for larger primes (`BIG_Q = 2⁶¹−1`, sound to `n ≲ 60`).
`--field {q,big,small}`. **See the "Step 6" section below.**

**Step 7 (S505): real leaf openings (`--pcs`) replace the per-layer `mle_eval`
stand-ins.** The `O(√x)` direct `mle_eval` leaf opens that grounded each carried
claim are replaced by sum-check OPENINGS (`leaf_open.py`: `open_eval`/`open_batch`)
whose residuals **thread to the chain base** instead of closing per layer. Threaded
through `line_batch_pair`/`batch_on_table` (both carried-claim folds) and the now-
redundant `s2`/`s_B` closes in `verify_affine_region`/`verify_trace_region`. The
chain verdict is unchanged (selftest case 20, over q & BIG_Q, automaton &
delegated+structured, composed with batch_trace+batch_wiring); `--bench-pcs` shows
a **stable ~1.2–1.36× verifier reduction**. This is a constant-factor win + a
working demonstration of the unconditional residual-threading architecture, NOT yet
asymptotic Õ(√x): the dominant per-layer `O(2^nb)` term that REMAINS is
`verify_trace_region`'s nb Ub-bit-table openings (line 429), whose residuals are
per-layer witness data needing the batched-trace integration. Default `pcs=False`
keeps prior artifacts verbatim. **See `leaf_open_results.md` for the full
treatment.**

**Step 8 (S506): the LAST per-layer `O(2^nb)` verifier term is gone (`--batch_ub`).**
After S505, exactly one per-layer `O(2^nb)` verifier operation remained:
`verify_trace_region`'s `nb` Ub-bit-table openings `Ub_{Lv-nb+k}~(r_C)` (the B2
routing check that pins `g1_trace` to the certified quotient `u_e`) — per-layer
*witness* data, so their residuals have no carried chain claim to thread to the base.
This cycle **defers** them (`ub_defer`, mirroring the `batch_wiring` accumulator: the
prover supplies the scalars, the verifier folds them into its `expect` in `O(nb)`) and
**discharges all `K·nb` of them in ONE degree-2 sum-check** against the SAME stacked Ub
cube `batch_trace`'s zero-test already commits. New primitive
`batched_trace.verify_ub_openings_batched`; the γ-RLC factorizes so the whole batch is
a single `claim = Σ_w B[w]·C[w]` over the `(Lk+nb)`-cube — `B` the verifier-anchored
per-layer eq weights, `C` the γ-fold of the committed low-`nb` Ub tables. Requires
`pcs` + `batch_trace`. **Measured asymptotic fact:** the verifier-side `O(2^nb)`
Ub-leaf eval count `ub_leaf_v` drops from **`K·nb` to exactly 0** at every `n` (moved to
the prover, `ub_leaf_p = K·nb`); the per-layer verifier is now **leaf-eval-free**, so the
only `O(√x)` verifier work left is one-time (the two `S₀` base closes + the batched
discharges). **The chain verifier is honestly Õ(√x) end-to-end** — closing the "Honest
scope" caveat that has stood since step 1. **Honest wall-clock scope:** the measured
`t_verifier` drop is only `~1.1–1.2×` and roughly flat in `n` (a vectorised numpy
`mle_eval` over an `O(√x)`-size array is cheap relative to the Python-loop wiring/eq
recomputes that dominate the measured `t_verifier`), so the `O(K·nb·2^nb)→`one-time
improvement does NOT surface as a growing wall-clock ratio at reachable `n`. This
**corrects the S505 NEXT-ACTION prediction** that the ratio would grow with `n`: it
grows in *op-count*, not in measured wall-clock until past the array-size crossover.
Verdict unchanged (selftest case 21, over q & BIG_Q; honest == sieve, delta_pi + liar
rejected). Default `batch_ub=False`. **See the "Step 8" section below.**

## Step 8 — batched Ub openings (S506): the per-layer verifier is leaf-eval-free

**The single remaining per-layer `O(2^nb)` verifier term (after S505).** `--pcs`
threaded every CARRIED-claim leaf open's residual to the `S₀` base. What it could not
touch: `verify_trace_region`'s B2 check evaluates `Ub_{Lv-nb+k}~(r_C)` for the `nb` low
bits of the certified quotient `u_e` (bit `nb−1−k` of `u_e`), `nb` direct `mle_eval`s
per layer = `O(K·nb·2^nb) = Õ(x)` over the chain. These are *witness* tables (no carried
claim to thread), so they need a commitment-shared discharge — and `batch_trace` already
stacks and commits exactly these Ub tables across all `K` layers in its zero-test.

**The construction (the γ-RLC factorization).** Each layer `l` emits the obligation
`(l, r_C^l, [ub_{l,k}]_{k<nb})`. With `γ ← F_q` and `β = γ^{nb}` the per-`(l,k)` weight
`β^l·γ^k = γ^{l·nb+k}` is a *distinct* power, so the weighting is a genuine RLC over all
`K·nb` openings. The combined identity

```
claim := Σ_{l,k} γ^{l·nb+k} ub_{l,k}  ==  Σ_w B[w]·C[w]
  B[w] = Σ_{l<K} β^l [w_L=l] eq~(r_C^l, w_e)     (verifier-anchored: public r_C, β)
  C[w] = Σ_{k<nb} γ^k Ub_{(nb−1−k)}[w]           (γ-fold of the committed Ub cube)
```

is ONE degree-2 sum-check over the stacked `(Lk+nb)`-cube (`Lk=⌈log₂K⌉`). The cross term
`Σ_w B[w]C[w] = Σ_{l,k} β^l γ^k Ub_{(nb−1−k)}^{(l)}~(r_C^l)` is exactly the γ-weighted
*true* openings, so an honest prover's `claim` matches; a single wrong `ub_{l,k}` flips
`claim` (caught by the sum-check round-1 identity, error `≤ (K·nb−1)/q + 2(Lk+nb)/q`).
**Verifier:** recompute `B~(r*)` in `O(K(Lk+nb))` (the anchor — `B` is public) and take
`C~(r*)` from the sum-check's folded scalar (the committed Ub-cube opening, exactly as
the zero-test trusts its U/R/Qv/bit scalars). One-time `O(K(Lk+nb)) = Õ(√x)`, replacing
the per-layer `O(K·nb·2^nb)`.

**Measured (`--bench-ub`, full config delegate+structured+batch_trace+batch_wiring+pcs,
q field):**

```
  n    K  nb   K·nb   vleaf_off  vleaf_on   tV_off    tV_on   tVr    tP_off    tP_on
  8    6   4     24        24        0      5.93m    5.57m  1.06x    8.29m     9.06m
 10   11   5     55        55        0     15.47m   13.45m  1.15x   17.67m    20.36m
 12   18   6    108       108        0     36.19m   31.02m  1.17x   35.31m    39.54m
 14   31   7    217       217        0     69.37m   56.66m  1.22x   63.72m    68.31m
 16   54   8    432       432        0    177.12m  157.27m  1.13x  146.38m   186.09m
```

`vleaf_on = 0` at every `n` is the headline: the per-layer verifier has no table-sized
eval. `tP` rises by the shifted openings. `comm` rises by one batched sum-check
(`1 + 3(Lk+nb)`, e.g. `+34` at `n=12`). The verdict and `claimed π(x)` are identical to
`batch_ub=False`.

## Step 9 — the whole-chain verifier op-count CURVE (S507): Θ(x) → Θ(√x)

Steps S505/S506 made the per-layer verifier **leaf-eval-free** (an op-COUNT fact at
a single `n`). What remained to turn the Õ(√x) end-to-end verifier claim into a
**measured curve**: a whole-chain field-op-count attribution that shows the *leading
term* fall from Θ(x) to Θ(√x). Wall-clock `t_verifier` cannot show this — the
dominant sum-check FOLDING is untimed prover work (S501), and vectorised numpy
flattens the asymptotics (S506, `tVr ≈ 1.1–1.2×` flat). The honest headline is an
op-count.

**The metric (`--bench-verifier-ops`).** A counter `_acct_vleaf` tallies every
VERIFIER **large-table evaluation** — a direct `mle_eval` of a committed `2^nb = √x`-size
table — weighted by table length, split into the **per-layer critical path** (scales
with `K = π(√x)`) vs the **one-time** terms (the two `S₀` base closes). This is the
*only* verifier operation whose per-call cost grows with `x`: every other verifier step
— sum-check round checks, `eq_point`/`le_eval`/`band_eval` point evals `O(nb)`, the
delegated wiring carry chain `O(nb·log p)` (S495), and the batched-discharge MLE
recomputes `O(K·polylog)` (S502/S503/S506) — is polylog or `O(K·polylog)`, hence bounded
by Õ(√x). So `vleaf_ops` **is** the leading term, and `comm` (printed) corroborates that
the entire non-leaf residual stays Õ(√x) in all configs.

Three configs, **all delegate+structured** (so the wiring verifier is already polylog and
the leaf opens are the sole `x`-scaling verifier term), differing only on the leaf-opening
axis:

```
  n        x    V  nb    K     pi |   a no-pcs: pl_cnt   total x4? |     b pcs: pl_cnt   total x4? | c pcs+bt+ub: pl_cnt total x4?
  8      255   15   4    6     54 |               53      880   -  |             24     416   -  |               0      32   -
 10     1023   31   5   11    172 |              109     3552 4.0x |             55    1824 4.4x |               0      64 2.0x
 12     4095   63   6   18    564 |              197    12736 3.6x |            108    7040 3.9x |               0     128 2.0x
 14    16383  127   7   31   1900 |              371    47744 3.7x |            217   28032 4.0x |               0     256 2.0x
 16    65535  255   8   54   6542 |              701   179968 3.8x |            432  111104 4.0x |               0     512 2.0x
 18   262143  511   9   97  23000 |             1357   695808 3.9x |            873  448000 4.0x |               0    1024 2.0x

fitted leading-term exponent  alpha  (total_ops ~ x^alpha; log2(total) vs n, last 4 pts):
        config  alpha_ops  alpha_comm
      a no-pcs      0.961       0.604
         b pcs      0.998       0.602
   c pcs+bt+ub      0.500       0.604
```

**The curve is the proof.** Per step (`Δn = 2`, i.e. `x` quadruples) configs (a)/(b)
grow `~4×` (∝ `x`), config (c) grows **exactly 2.0×** (∝ `√x`). The fitted exponents:
(a) `α_ops = 0.961`, (b) `0.998` — `Θ(x)`; (c) `α_ops = 0.500` — `Θ(√x)`. The exact
per-config leaf-eval COUNTS (asserted in the bench and selftest 22) are
`K·(nb+5)−1` (a), `K·nb` (b, the Ub openings only), and **0** (c); the one-time term is
`2·2^nb` in all three. Config (c)'s total `vleaf_ops` is *only* the two `S₀` base
closes — `2·2^nb` — so the entire Θ(x) per-layer term that (a)/(b) carry has collapsed
to a one-time Θ(√x) term.

**Honest scope.** (i) The metric counts large-table (`2^nb`) evals only — the leading
term. Config (c)'s *total* verifier op-count is `vleaf_ops + non-leaf work`; the non-leaf
work (sum-check checks + batched-discharge recomputes) is `O(K·polylog) = Õ(√x)`, proven
by construction in S495/S502/S503/S506 and tested there, and reflected here in the
`comm` column (slope `α_comm ≈ 0.60` — `√x` with polylog factors above it, sub-linear in
all three configs, so it never re-introduces a Θ(x) term). So config (c)'s grand total is
`Õ(√x)` (effective exponent ~0.6 at these `n`, polylog above `√x`), not literally `x^{0.5}`
— which is exactly the "leading term `Θ(√x·polylog)`" the milestone claims. (ii) `α_ops`
for (a) reads `0.961` not `1.000` because `total(a) = (K(nb+5)+1)·2^nb` carries the
`π(√x)/√x` and `nb` log-factors that bend the finite-`n` slope slightly below 1; the
per-step `~4×` and (b)'s `0.998` confirm the term is genuinely linear. (iii) Field
(`q` vs `BIG_Q`) does not change the op COUNT — same tables, same number of evals
(selftest 22 asserts this); the bench uses default `q` for speed.

The remaining Θ(√x) verifier work is now entirely **one-time** (the two `S₀` base closes,
`2·2^nb`, plus the `O(K·polylog)` batched discharges) — so a real outer multilinear
commitment on those two base closes would make the whole verifier *succinct* (polylog)
under a hashing assumption. The unconditional verifier is otherwise complete at Õ(√x).

## Step 6 — the field lift through the compressed chain (S499)

The hole closed by S498 lived in `compressed_prover_mult_trace.py`'s isolated
primitive; the chain (`compressed_layer.py`) and the delegated wiring
(`lucy_dp_delegated_wiring.py`) still hard-coded the global `Q = 2³¹−1` imported
from `lucy_dp_verification` (every `% Q`, `pow(.,Q−2,Q)`, `eq_table`, `sumcheck`,
and each automaton eval). S499 threads a `q=` parameter through all of them,
exactly as S498 did for the primitive:

- **Helpers reused, not duplicated.** Both files now import the
  field-parameterised `lagrange_eval / eq_table / eq_point / mle_eval / sumcheck`
  (plus `_dt / _asum`) from `compressed_prover_mult_trace`, which default to
  `q = DEFAULT_Q`. The two pure-Python automata `ge_const_eval` / `w_div_eval`
  in `lucy_dp_verification` gained a `q=Q` kwarg. `waff_eval`, `le_eval`,
  `band_eval`, and the whole reduction/wiring stack carry `q` down.
- **Dtype.** `_dt(q)` picks `uint64` for `q ≤ 2³¹−1` (a+b of two products stays
  below 2⁶⁴) and `object` (exact Python ints) above it; `_asum` reduces either.
  Index/witness arithmetic stays `int64` (safe for `x < 2⁶³`, i.e. `n ≲ 62`).
- **Backward-compatible.** Default `q = Q` → the entire transcript is verbatim:
  all prior selftests (incl. the structured-vs-dense bit-identity, selftest 15)
  pass unchanged.

**Verified (selftest case 16):**
- The whole compressed chain runs over `BIG_Q` (object dtype) at `n ∈ {8,10,12}`,
  all three modes (automaton / delegated-dense / delegated-structured):
  `claimed π(x) == sieve π(x)` and every chain cheat (`delta_pi`,
  self-consistent liar) rejected. This proves the q-threading is arithmetically
  correct end-to-end over a larger-**characteristic** prime.
- The lift **closes the wrap-around aliasing hole in the chain's exact
  trace-region config** `build_witness(x, p, nb, dstart=1)`: at `n = 22`
  (`x = 2²²−1 > SMALL_Q = 1000003`) a forged quotient `u''` aliases `X mod
  SMALL_Q` (so the too-small prime **accepts** it, 1/1) but violates
  `U·a + R − X = 0` over `BIG_Q` (which **rejects** it, 5/5). Same `forge_alias`
  / `verify_constraints` machinery S498 validated, now exercised in the chain's
  configuration.

**Honest scope — what n>31 means here.** The *soundness* failure of the demo
prime at `n>31` lives entirely in the trace certifier (`verify_constraints`),
which is field-parameterised and tested directly above at the chain's
`build_witness` config. A *full* `run_chain` at `n>31` is **prover-bound, not
field-bound**: `compressed_lucy` is an `O(x/log x)` pure-Python DP (≈10⁹
iterations at `n=34`), and the per-layer cubes are `2^{bits(√x)} ≈ √x` object
arrays over `K ≈ π(√x) ≈ 12000` layers — minutes-to-hours regardless of the
prime. So the n>31 evidence is at the trace certifier (the seat of the hole),
while whole-chain correctness over the lifted field is shown at small `n`. A
numpy Mersenne `mulmod` for `2⁶¹−1` (32-bit-limb schoolbook + the `mod 2⁶¹−1`
fold) would restore uint64-class speed for the object path; it is a pure
speed change, no structural difference, and is the next step if a large-`n`
chain benchmark is wanted.



Step 4 met the verifier target (`Õ(√x)`) but, by delegating, moved the binding
cost to the **chain prover**: each delegated wiring's backward sweep built five
`2^{2l}`-size tables (`inner_r_table` + tiled `EQ/LT/V`) and ran a `2l`-variable
sum-check over them — `O(nb·p²)` per Lucy layer, `Σ_{p≤√x} nb·p² ≈ Õ(x^{3/2})`.
S496 showed this `p²` is an *implementation* artifact: the summand
`F(σ',σ)=E(σ')P(σ')R_j(σ',σ)L(σ)v_j(σ)` factors at the `σ'/σ` boundary, `R_j` is
the sparse width-4 affine operator (4 nonzeros/column) and `v_j` has support
`≤p`, so splitting the sum-check into two `O(p)`-size phases gives an **identical
transcript** in `O(nb·p)` prover work (`chain_layer_reduce_structured`).

This cycle wires that into the compressed chain. The `structured` param is plumbed
through `verify_waff_value → verify_affine_region → large_reduce` and through
`small_reduce`, reaching the existing S496-ready
`inner_verify_div(..., structured=...)`. Both wirings the compressed chain
delegates are covered: the **small division** (`accept_rem=None`) and the **large
affine image** (`accept_rem=p−1`, the `verify_waff_value` inner claim). The
comparator fold in `verify_waff_value` is already `O(2^nb)` (p-independent) and is
untouched; only the delegated chain sweep changes.

**Drop-in is exact (selftest case 15).** Because the transcript is bit-identical,
`structured=True` vs `False` must agree on EVERY observable: selftest 15 asserts
identical accept/reject + comm for the two delegated primitives, identical
accept/reject + comm for the full two-sided layer across honest + 10 cheat
classes, and identical accept/reject + comm + claimed π(x) for the full chain
(honest, `delta_pi`, self-consistent liar) at `n∈{8,10,12}`. Default stays
`structured=False` so all existing artifacts and the dense path are unchanged.

**Measured — the wiring prover flattens `p²→p` (`--prover-bench`, fixed nb=14,
per-evaluation wall-clock ms, dense vs structured, both wirings):**

| p | l | div dense | div struct | ratio | dense/p² | struct/p | aff dense | aff struct | ratio |
|---|---|-----------|-----------|-------|----------|----------|-----------|-----------|-------|
| 7    | 3  | 22.9   | 18.3 | 1.25  | 467.0 | 2607 | 23.5   | 18.5 | 1.27  |
| 13   | 4  | 31.6   | 24.7 | 1.28  | 186.8 | 1897 | 31.2   | 24.5 | 1.27  |
| 31   | 5  | 44.2   | 29.6 | 1.50  | 46.0  | 953  | 43.0   | 29.4 | 1.46  |
| 61   | 6  | 60.9   | 35.9 | 1.70  | 16.4  | 588  | 59.7   | 35.2 | 1.70  |
| 127  | 7  | 107.6  | 40.6 | 2.65  | 6.67  | 319  | 90.2   | 34.2 | 2.63  |
| 251  | 8  | 226.0  | 40.2 | 5.62  | 3.59  | 160  | 234.6  | 51.4 | 4.57  |
| 509  | 9  | 951.6  | 51.7 | 18.42 | 3.67  | 102  | 774.1  | 62.8 | 12.33 |
| 1021 | 10 | 3420.1 | 74.0 | 46.20 | 3.28  | 73   | 3353.1 | 59.8 | 56.12 |

The **`dense/p²` column flattens to ~3.3** (dense is `Θ(p²)`); the structured
column grows far slower (dense ×149 from p=7→1021, structured ×4), so the **ratio
climbs 1.25× → 46×** — bit-for-bit mirroring S496's standalone `1.3×→43.9×`. The
`struct/p` column *decreases* (the structured per-layer work is `O(p)` but
overhead-dominated at the demo's tiny `l`, exactly as S496 noted), confirming the
structured prover is sub-quadratic — `O(nb·p)` is an upper bound, not yet the
binding term at these sizes. **Whole-chain wiring prover: `Σ_{p≤√x} nb·p² ≈
Õ(x^{3/2}) → Σ nb·p ≈ Õ(x)`** — the full compressed delegated π(x) chain prover
is now Õ(x), matching the substrate prover, with the verifier still Õ(√x).

**Open item 1 is complete end-to-end on the REAL compressed chain:** verifier
Õ(√x) (delegated, no p-linear term) + prover Õ(x) (compressed √x-state substrate
+ structured wiring sweep), unconditional, both wirings delegated. The remaining
work is the secondary list: field lift `|F|>x`, a real outer poly-commitment, and
the large-x benchmark.

## Step 4 — wiring delegation (S495): removing the last p-linear verifier term

Open item 1's verifier had ONE super-polylog term left: per layer the verifier
ran a long-division automaton with a running remainder in `[0,p)` to evaluate
the wiring MLE — `w_div_eval` (small side) and `waff_eval` (large affine side),
both `O(nb·p)`. Summed over all primes that is `Σ_{p≤√x} nb·p ≈ nb·x/ln x =
Õ(x)` — the verifier was secretly **linear in x**, not √x. Step 4 delegates
both wirings to the prover via the carry-chain GKR template already built in
`lucy_dp_delegated_wiring.py`, bringing the per-layer wiring check to
`O(nb·log p)` and the whole-chain verifier to `Σ_{p≤√x} nb·log p ≈
nb·θ(√x)/ln2 = Õ(√x)` (Chebyshev: `Σ_{p≤√x} log p = θ(√x) ~ √x`).

**Small side (drop-in, as the note predicted).** `small_reduce`'s phase-B
wiring `W(v,u)=[u=⌊v/p⌋]` is *exactly* the division relation
`lucy_dp_delegated_wiring.inner_verify_W` already delegates (affine relation
`s'=2s+a−bp`, width-4 carry automaton `AFF~_c`, a `2l`-variable sum-check per
bit, `l=⌈log₂p⌉`). The prover supplies the scalar `W~(r_v,r_u)` (computed in
Õ(√x) with **no** automaton — a single `Σ_v eq(r_v,v)eq(r_u,⌊v/p⌋)` cube
contraction), the verifier binds it to the phase-B `finalB` and proves it
through the chain.

**Large side (the genuinely new piece).** `waff_eval` is the affine-image
relation `[e'=p·e+(p−1)]` **crossed with a comparator `[e≤ecut]`**, so it is
not a pure division. The MLE of a product `C(e)·D(e,e')` is *not* the product
of MLEs, so the comparator cannot simply be factored out. Decomposition used
(`verify_waff_value`): group by `e`,
`waff_eval(r_v,r_u) = Σ_e eq(r_v,e)·C(e)·Φ(e)`, `C(e)=[e≤ecut]`,
`Φ(e)=eq(r_u,(e+1)p−1)`. An **inner sum-check over the e-cube** (degree 3)
folds the comparator **as-is** — its MLE is the existing `O(nb)` `le_eval`, no
new automaton — and the eq-factor, leaving one claim `Φ~(r_e)=D~(r_e,r_u)` for
the **pure** affine-image relation `D(e,e')=[e'=p·e+(p−1)]`. That is delegated
to the SAME carry chain with the dividend `=e'` (point `r_u`), the quotient
`=e` (point `r_e`), and **acceptance restricted to final remainder `p−1`**
(`inner_verify_div(accept_rem=p−1)`) instead of the division wiring's "any
remainder". The only change to the carry-chain template is this accept
weighting in its sum-check #0 (`[σ=p−1]`, MLE `eq_point`, vs `[σ<p]`, MLE
`lt_const_eval`); the backward chain sweep — whose `LT~` masks junk states in
the *agreed extension*, unrelated to acceptance — is byte-for-byte identical.

**Measured — the wiring check flattens in p (`--wiring-bench`, fixed nb=12,
random points, verifier-time ms):**

| p | l=⌈log₂p⌉ | waff automaton | waff delegated | div automaton | div delegated |
|---|---|---|---|---|---|
| 3 | 2 | 0.101 | 0.877 | 0.048 | 0.787 |
| 13 | 4 | 0.293 | 1.417 | 0.127 | 1.390 |
| 61 | 6 | 0.848 | 2.045 | 0.398 | 2.060 |
| 127 | 7 | 1.790 | 2.118 | 0.691 | 2.156 |
| 251 | 8 | 3.178 | 2.764 | 1.171 | 2.770 |
| 509 | 9 | 4.387 | 2.948 | 2.053 | 2.747 |
| 1021 | 10 | 6.258 | 3.102 | 2.839 | 3.248 |

The **automaton columns grow ~linearly in p** (waff ×62 from p=3→1021); the
**delegated columns track l=⌈log₂p⌉** (waff ×3.5 over the same range, i.e.
~log p). waff crosses over at **p≈251**, the division wiring near **p≈1021**;
beyond that the delegated verifier is unboundedly cheaper (at `p≈√(10¹⁰⁰)`
the automaton's `p` is astronomical, the chain's `log p` is ~166 bits).

**Honest trade-offs (both real, both reported):**
- **Verifier wins asymptotically, loses at demo sizes.** The chain's
  per-layer constant (several `2l`-variable sum-checks) is larger than one
  automaton pass, so the *full chain* delegated verifier is ~2.3× **slower**
  in wall-clock at `n≤16` (`--bench`: e.g. n=16 tV 88→212 ms), because every
  prime there is small (`p≤V≤255`, mostly below the crossover). The deliverable
  is the **scaling** — `O(nb·log p)` per layer, no p-linear term — verified by
  `--wiring-bench`; at the demo's tiny x, Õ(x) and Õ(√x) are empirically
  indistinguishable (both grow ~2.2–2.4× per √x-doubling).
- **Communication rises ~5×** (n=16: 15 208 → 87 226 elems): the chain
  transcripts. Still `Õ(√x)`-compatible (`O(K·nb·log p)` field elements).
- **Prover cost rises** — and this is the **remaining bottleneck**. The current
  `lucy_dp_delegated_wiring` chain prover materializes dense `2^{2l}`-size
  tables per chain layer (`inner_r_table`, the tiled `EQ/LT/V`), i.e.
  `O(nb·p²)` per Lucy layer, `Σ_{p≤√x} nb·p² ≈ Õ(x^{3/2})`. So delegating the
  wiring trades a verifier `Õ(x)→Õ(√x)` for a prover `Õ(x)→Õ(x^{3/2})`. The
  `O(p²)` is an *implementation* artifact, not fundamental: `R~_j` is the
  sparse width-4 affine operator (4 nonzeros per state) and `v_j` has support
  `≤p`, so a structured prover that never densifies is `O(nb·p·l)` per layer
  (`Õ(x)` total). Wiring that structured prover is the clean next step; the
  **verifier target of open item 1 is met now**.

**Soundness — two NEW delegation-specific liars per side, all rejected
(12/12 at n=14, 6/6 across selftest):**

| cheat | mechanism caught |
|---|---|
| `waff_forge` / `small_forge` | corrupt the wiring table **and** forge the matching scalar `=finalB·s2⁻¹` (so the algebra check passes) → the inner sum-check's round-1 sum over the **honest** Φ/W tables ≠ the forged claim |
| `waff_chain` / `small_chain` | honest scalar, but the prover corrupts the delegated chain mid-way (consistent re-propagation) → the chain's base check `v~_0 = e_0` fails |

The forge liar is the load-bearing one: it proves you cannot make the verifier
accept a wrong wiring **even when you control the claimed scalar** — the value
must be re-derived by the verifier-checkable inner protocol. All 9 step-3
cheats (`wrong_C/G/split/u/waff/small_*/carry/delta_pi/self-consistent-liar`)
remain rejected in `--delegate` mode (13 layer + chain classes total).

**Cross-check (selftest case 13):** the delegated protocols prove the **same
scalar** the automata compute — `verify_waff_value(waff_eval(r_v,r_u))` and
`inner_verify_div(w_div_eval(r_v,r_u))` accept, `…+1` and a lying chain reject,
at random field points for several `(nb,p,ecut)`.

## Step 3 — small-side layer + chaining (S494)

**Small-side layer** (`small_reduce`). A claim `C=S̃_i^small(z)` (value-addressed
`v∈[1,V]`) reduces by the **same two-phase shape** as the base protocol's
`verify_layer`, but over the `2^nb≈√x` compressed cube:
- **phase A** (deg ≤ 3): `C = Σ_v eq(z,v)[S_{i-1}^small(v) − B(v)(G(v)−(i−1))]`,
  `G(v)=S_{i-1}^small(⌊v/p⌋)`. The wiring `v→⌊v/p⌋` is **value→value, entirely
  inside the small side** — no trace, no cross-side lookup.
- **phase B** (deg 2): `g1 = Σ_u W(r_v,u)·S_{i-1}^small(u)`, `W(v,u)=[u=⌊v/p⌋]`
  verified by the existing **long-division automaton** `w_div_eval` (O(nb·p)).
- emits **two** point-claims on `S_{i-1}^small`: `s1@r_v` and `s2@r_u`.

**The gate is a BAND, not a bare comparator — a correction to the planned
design.** The NEXT-ACTION note said `B(v)=[v≥p²]` is "a comparator, no
validity-fold needed (small keys never leave [1,V])". That holds **only when
`V=2^nb−1`** (i.e. `x=2^n−1`, `n` even → no padding). For general `V` the cube
has padding columns `(V,2^nb)` that are 0 every layer; with a bare `[v≥p²]`
those columns would receive a nonzero correction `−(S(⌊v/p⌋)−(i−1))≠0` and the
prover's update MLE would no longer equal the true padded `S_i^small`. The
correct gate is the **band** `B(v)=[p²≤v≤V]` (`band_eval`); when `V=2^nb−1` it
collapses to the comparator, so the note is the no-padding special case. A
second boundary surfaced: when `p²>V` the band is **empty** and `band_eval`
must return 0 (the `le(hi)−le(lo−1)` form is wrong for `lo>hi`) — now guarded
and exhaustively selftested.

**The two-sided layer & chain** (`run_two_sided_layer`, `run_chain`). We carry
a **width-2 claim vector** `(large-claim about S_i^large, small-claim about
S_i^small)`. Per layer:
- `large_reduce` (step 2) consumes the large-claim → a large-claim about
  `S_{i-1}^large` (line-batched `s1@r_v`+`s2_aff@r_u`) **and** a small-claim
  `s_B@r_B` about `S_{i-1}^small` (its trace region);
- `small_reduce` (skipped only at the top layer `i=K`, where no `S_K^small` is
  referenced) consumes the small-claim → two small-claims about `S_{i-1}^small`;
- the **same-table** small claims (`s_B` + the two from `small_reduce`) fold to
  ONE by sequential **line restriction** (`batch_on_table`); cross-table (large
  vs small) claims stay separate — exactly the "carry a small fixed-width claim
  vector" plan.
The chain starts from a single Boolean claim `S_K^large(e=0)=π(x)` and runs to
`S_0`, opened on both sides (`mle_eval`, the leaf stand-in — `O(√x)` once,
within budget; the `S_0` bases have no cheap closed form: `⌊x/(e+1)⌋` large,
padded `max(v−1,0)` small).

**Measured (`q=2³¹−1`, `run_chain`, full K-layer chain):** `claimed π(x) ==
sieve π(x)` at every `n`.

| n | V=√x | nb | K=π(√x) | π(x) | t_prover (ms) | t_verifier (ms) | comm (elems) |
|---|------|----|---------|------|---------------|-----------------|--------------|
| 8  | 15  | 4 | 6  | 54   | 4.3   | 3.6  | 724   |
| 10 | 31  | 5 | 11 | 172  | 10.0  | 8.4  | 1 746 |
| 12 | 63  | 6 | 18 | 564  | 22.3  | 17.8 | 3 562 |
| 14 | 127 | 7 | 31 | 1900 | 48.1  | 38.5 | 7 406 |
| 16 | 255 | 8 | 54 | 6542 | 109.6 | 93.2 | 15 208 |

Both prover and verifier grow ≈2.2× per `+2` in `n` (one √x-doubling) at these
sizes; comm ≈ `K·O(nb²)`. **Asymptotics (structural, not an empirical fit of 5
points):** per layer the prover works only the `2^nb≈√x` compressed cube (it
never materializes the `2^n≈x` table), so the **prover is Õ(x) total** =
`K≈√x/ln√x` layers × Õ(√x) per layer — the open-item-1 target (vs the demo's
Õ(x^{3/2})). The **verifier is Õ(√x) per the base-protocol accounting** once
`waff`'s O(nb·p) is delegated to O(nb·log p) (the only super-polylog term;
`lucy_dp_delegated_wiring`'s carry chain already does this for the analogous
division wiring — deferred, not a new obstruction).

**Soundness (all cheats rejected; 20/20 at n=12, 8/8 across selftest):**

| cheat | mechanism caught |
|---|---|
| `delta_pi` (claim π(x)+1) | layer-K phase-A round-0 sum (honest `S_{K-1}` tables) |
| self-consistent liar @ i₀ (corrupt one DP entry, re-propagate up) | layer i₀: corrupted carried claim ≠ phase-A sum over the honest `S_{i₀-1}` below |
| `small_C` / `small_G` / `small_u` | small phase-A sum / phase-B division anchor (`g1≠ΣW·S`) / `w_div_eval` mismatch |
| `carry` (tamper the cross-side `s_B` before batching) | `batch_on_table` endpoint check `h(0)=mle_eval(S,r_B)≠s_B` |
| `wrong_C`/`G`/`split`/`u`/`waff` (large side) | as in step 2 (unchanged) |

The self-consistent liar is the load-bearing test: it confirms the **cross-layer**
binding works — because the small side's recurrence never reads the large side,
corrupting a large entry leaves all small layers honest, so the catch lands on
the large side at exactly layer i₀ (verified at i₀ ∈ {1, K/2, K}).

## What is built

1. **The compressed Lucy_Hedgehog DP** (`compressed_lucy`) — the prover
   substrate that was missing (the existing protocols materialize the full
   `2^n≈x` table, an Õ(x^{3/2}) demo prover). State, `V=⌊√x⌋`:
   `small[v]=S(v)` value-addressed `v∈[1,V]`; `large[e]=S(⌊x/(e+1)⌋)`
   d-addressed `e=d-1∈[0,V-1]`. Validated: `π(x)=large[0]` equals a direct
   sieve for every `n∈[6,20]`.

2. **One large-side compressed layer verifier** (`run_layer`). A claim
   `C=S̃_i^large(z)` reduces by
   - **phase A** (deg ≤ 3): `C = Σ_e eq(z,e)[S_{i-1}^large(e) − B(e)(G(e)−(i−1))]`
     with `B(e)=[e≤Bcut]` a **single threshold** that folds validity `e≤V−1`
     **and** the `p²` gate (`⌊x/(e+1)⌋≥p² ⇔ e≤⌊x/p²⌋−1`); verifier evaluates
     it with one `le_eval`. Output claims `s1=S̃_{i-1}^large(r_v)`, `g1=G̃(r_v)`.
   - **phase B**: `g1 = g1_aff + g1_trace` (prover sends the split, verifier
     checks the sum), `ecut=⌊V/p⌋−1`:
     - `g1_aff = Σ_{e≤ecut} eq(r_v,e)·S_{i-1}^large((e+1)p−1)` — the
       **affine index map** `e↦p·e+(p−1)` (target `⌊x/((e+1)p)⌋` is itself a
       large key when `(e+1)p≤V`), verified by `waff_eval`: a long-division
       automaton on `e'` (quotient bits forced to equal `e`, accept on final
       remainder `p−1`) crossed with a `e≤ecut` comparator. **No trace.**
     - `g1_trace = Σ_{ecut<e≤V−1} eq(r_v,e)·S_{i-1}^small(⌊x/((e+1)p)⌋)` — the
       **step-1 trace primitive**: `build_witness`+`verify_constraints` certify
       every `u_e=⌊x/((e+1)p)⌋` (degree-3 zero-test, no division), then a
       **band-masked** B1/B2 routes the small-side lookup.
   - **line batching**: `s1@r_v` and `s2_aff@r_u` (both on `S_{i-1}^large`)
     collapse to one claim (the existing line-restriction trick); the trace
     produces `s_B=S̃_{i-1}^small(r_B)`, carried for the small-side layer.

## Two design facts that keep the sizes at √x

- **The band mask is load-bearing and verified.** The trace constraint test
  certifies `u_e` over the *whole* e-cube (so the verifier's wiring
  `a~(e)=(e+1)p` stays a clean degree-1 line and `e=0` never divides by zero).
  The affine/padding columns have `u_e>V`; the band mask `T(e)=[ecut<e≤V−1]`
  zeroes them out of the routing, and its MLE is checked in B2 by
  `band_eval = le(hi) − le(lo−1)` (the MLE is linear over that pointwise
  Boolean difference — no extra automaton). Without the verified mask a prover
  could route an affine column's low bits into a real small bucket.
- **The lookup is over nb bits, not Lv bits.** `build_witness` needs
  `Lv≈bits(⌊x/p⌋)` bits to *certify* the large affine `u_e`; but the **routing**
  takes the **low `nb=bits(V)` bits** of the certified `u_e` (trace `u_e≤V` so
  its high bits are 0), so `omega`/`S_small` stay size `2^nb≈√x`. The
  constraint test itself runs on the e-cube (`2^nb≈√x` rows) with `Lv=O(log x)`
  bit tables → prover Õ(√x).

## Measured (q = 2³¹−1)

Honest layer accepted; **output point-claims verified equal to the true MLEs**
of `S_{i-1}^small` (at `r_B`) and `S_{i-1}^large` (the line-batched point) —
the reduction is consistent, not just internally accepting. All cheat classes
rejected (20/20 at n=12, and 8/8 across the selftest's primes/layers):

| cheat | mechanism caught |
|---|---|
| `wrong_C` (claimed `S_i^large(z)` +1) | phase-A round-0 sum |
| `G` (corrupted wiring table) | `g1 ≠ g1_aff+g1_trace` |
| `split` (g1_aff+1, g1_trace−1, sum preserved) | region sub-protocol claim |
| `u` (wrong trace quotient) | step-1 multiply+remainder zero-test |
| `waff` (wrong affine image index) | `waff_eval` mismatch |

**Scaling** (`--bench`, middle layer; `x` grows 256×, n 10→18):

| n | V | nb | layer | p | t_prover (ms) | t_verifier (ms) | comm |
|---|---|----|-------|---|---------------|-----------------|------|
| 10 | 31 | 5 | 5 | 11 | 0.53 | 0.621 | 116 |
| 12 | 63 | 6 | 9 | 23 | 0.66 | 0.680 | 145 |
| 14 | 127 | 7 | 15 | 47 | 0.74 | 0.902 | 176 |
| 16 | 255 | 8 | 27 | 103 | 0.89 | 1.252 | 209 |
| 18 | 511 | 9 | 48 | 223 | 1.16 | 2.629 | 244 |

`comm = O(nb²) = O(log²x)`; prover ≈ flat (works the `2^nb≈√x` cube). The
verifier's mild super-polylog growth is **entirely `waff_eval`'s O(nb·p)**
automaton (`p≈223` at n=18) — the *base*-protocol wiring cost, which
`lucy_dp_delegated_wiring` already showed reduces to O(nb·log p) by the
carry-chain. Wiring delegation of `waff` is deferred (next step), not a new
obstruction.

## Honest scope — what is and isn't done

- **Done:** the compressed DP substrate; one large-side layer's full
  reduction with the small/large dispatch; the affine wiring automaton; the
  band-masked trace lookup reusing step 1; line batching of the two
  large-side claims; output claims checked against true MLEs; 5 cheat classes.
- **Leaf closure** (`S` at `r_B`, `Ub_j` at `r_C`, the line-batched large
  point) is direct MLE evaluation — the same harness stand-in as step 1 for a
  poly-commitment opening / outer line-batching. The two `Ub_j` point-claims
  (`r_A` from the constraint test, `r_C` from B2) reconcile by the same trick.
- **Field caveat** (inherited from step 1): demo `q=2³¹−1` works because
  `u_e·a_e ≤ x < q` at the tested `n`; production lifts `|F|>x` or uses a
  small-field carry trace. Unchanged here.
- **DONE in step 3 (S494):** the **small-side** layer (value→value division
  wiring); **chaining** the full `K=π(√x)`-layer two-sided chain from a single
  `S_K^large(e=0)=π(x)` claim to `S_0`, `claimed π(x)==sieve π(x)`;
  self-consistent-liar soundness.
- **DONE in step 4 (S495):** **delegating both wirings** (`waff` and the small
  side's `w_div_eval`) from O(nb·p) to O(nb·log p) via the
  `lucy_dp_delegated_wiring` carry chain — the last super-polylog verifier term
  is gone, the verifier is cleanly **Õ(√x)** (`--delegate`). See "Step 4".
- **DONE in step 5 (S496→S497):** the **structured chain prover** — the
  delegated chain's dense `2^{2l}` tables replaced by the sparse width-4 affine
  operator so the *prover* drops from `Õ(x^{3/2})` back to `Õ(x)`. See "Step 5".
- **DONE in step 6 (S499):** the **field lift** — a prime-modulus parameter `q`
  threaded through the whole chain and the delegated wiring, defaulting to `Q`
  (verbatim) with `BIG_Q = 2⁶¹−1` (object dtype) removing the `u·a<q` shortcut;
  the chain's trace certifier now rejects the wrap-around alias the demo prime
  admits past `n≈30`. See "Step 6".
- **DONE in step 7 (S505):** the carried-claim leaf opens (the line-batched large
  point, the `s2`/`s_B` reconciliations) are real sum-check **openings** whose
  residuals thread to the `S₀` base (`--pcs`), unconditionally. See `leaf_open_results.md`.
- **DONE in step 8 (S506):** the **`Ub_j` openings** — the only per-layer `O(2^nb)`
  verifier term `--pcs` could not thread (per-layer witness data) — are deferred and
  **discharged in ONE sum-check against the committed stacked Ub cube** (`--batch_ub`,
  requires `pcs`+`batch_trace`). Verifier-side `O(2^nb)` Ub-leaf eval count `K·nb → 0`;
  the per-layer verifier is **leaf-eval-free**, the end-to-end verifier honestly **Õ(√x)**
  (the only `O(√x)` left is the two one-time `S₀` base closes). See "Step 8" above.
- **NOT done (next steps):** (c) a **real outer poly-commitment** to replace the
  two one-time `S₀` base `mle_eval` closes (computational, not unconditional); (d) the end-to-end
  **large-x benchmark** vs a real Lucy recomputation — see S500 below: the
  numpy Mersenne `mulmod` is built and threaded through the chain, but it does
  **not** unblock a large-`n` chain run; the chain wall is the *count* of small
  per-layer field ops (not the Lucy DP, and not the per-multiply cost).

## Fast BIG_Q path + where the chain time actually goes (S500)

The S499 NEXT ACTION assumed a large-`n` `run_chain` is bound by "object-dtype
arithmetic **and** the `O(x/log x)` pure-Python `compressed_lucy` DP". Measuring
both **corrects this**:

- **The Lucy DP is negligible.** `compressed_lucy(x)` is `0.5 ms` (n=16),
  `232 ms` (n=24), `868 ms` (n=26) — **< 0.1 % of `run_chain`**, which is
  essentially 100 % field arithmetic (n=18 chain: DP `4.6 ms` vs chain `48 s`).
  So speeding the DP is **not** the lever; the field work is.
- **The Mersenne `mulmod` is built and threaded chain-wide, bit-identically.**
  `_mul61`/`_sum61` (in `compressed_prover_mult_trace.py`) reduce
  `(a·b) mod (2⁶¹−1)` on uint64 via the `2⁶¹≡1` fold; the chain's array products
  route through `fmul` / `_modmul_arr` and the `np.add.at` accumulations
  (`build_omega`, `large_reduce`, `small_reduce`) fall back to exact object
  accumulation under the fast flag. Selftest §17 cross-checks the **whole big
  chain fast-vs-object, bit-for-bit** — same claimed π, accept/reject, comm, in
  all three modes, with chain cheats and the alias still rejected.
- **But the fast path is SLOWER for the chain** (`--bench-big`, delegate +
  structured, BIG_Q, uint64 vs object):

  | n | V=√x | K=π(√x) | obj ms | fast ms | speedup |
  |---|---|---|---|---|---|
  | 8 | 15 | 6 | 326 | 1476 | 0.2× |
  | 12 | 63 | 18 | 2134 | 9106 | 0.2× |
  | 16 | 255 | 54 | 15473 | 51404 | 0.3× |

  **Mechanism:** `_mul61` is ~24 vectorised numpy ops per multiply, so it only
  amortises on **large** arrays (≳1000 elems; the certifier crosses over at
  `D≈2000`, see `compressed_prover_mult_trace_results.md`). The chain is the
  opposite shape — `π(√x)` layers each doing **many small** cubes (nb-bit and
  `p`-size, ≤√x and usually ≪√x) — so at reachable `n` (√x ≤ 255) it is
  *operation-count-bound on small arrays* and the object path's 2-op multiply
  wins. The chain's own crossover is at large `n`, where the per-layer
  `√x`-cube sum-checks dominate, but that is beyond cheap object-reference
  measurement (n=24 object chain ≈ 20 min extrapolated).
- **Consequence for the benchmark:** the right lever is **reducing the number
  of small field-ops** (cross-layer / cross-cube batching so the fast path's
  width is exercised), not the per-multiply cost and not the DP. The mulmod is
  kept as an **opt-in** (`--fast-big`); BIG_Q stays on the object reference by
  default because it is faster at reachable `n`. The clean, measured win of the
  mulmod is the **certifier** (one large `√x`-cube: 1.7–3.8× at n=24–28).

## Cross-cube batching: trace (S502) + wiring (S504) → fast path wins (S504)

The S500 mechanism named the lever: **reduce the number of small field-ops**
(widen them by cross-cube batching) so the fast Mersenne path is exercised. Two
batched primitives were built standalone — the trace zero-test
(`batched_trace.py`, S502, wired via `run_chain(batch_trace=True)`) and the
wiring delegation (`batched_wiring.py`, S503) — and this cycle (**S504**) wired
the wiring batch end-to-end via `run_chain(batch_wiring=True)`.

**Defer-and-batch.** With `batch_wiring=True` (delegate only), the per-layer
delegated wiring chains are no longer discharged inline: a `defer` accumulator
threaded through `large_reduce`/`verify_affine_region`/`verify_waff_value`
(large affine image, `accept_rem=p−1`) and `small_reduce` (small division,
`accept_rem=None`) **collects** each obligation `(p, r_v, r_u, accept_rem,
claim, lie)` — the claim being the scalar the layer's outer sum-check already
pinned — and all `2K−1` obligations discharge in ONE
`batched_wiring.verify_wiring_obligations` after the layer loop. The
`verify_waff_value` `[e≤ecut]` comparator fold (O(2^nb), p-independent) stays
per-layer. Default `batch_wiring=False` keeps every prior artifact verbatim.

**Verdict unchanged** (selftest §19/§19b, over `q` & `BIG_Q`, structured &
dense, alone and composed with `batch_trace`): honest accepts & matches the
sieve; `delta_pi` and the self-consistent liar rejected; the wiring-specific
liars (`small_forge` at sum-check #0; `small_chain`/`waff_chain` in the batched
backward sweep; `waff_forge` in the per-layer comparator fold) rejected through
the batched discharge. Transcript is **not** bit-identical (the batch draws rng
after the layer loop), so only verdict/claimed are asserted, not comm.

**End-to-end headline** (`--n 16 --bench-combined`, BIG_Q, delegate+structured):

  | config | wall (ms) | comm | vs baseline |
  |---|---|---|---|
  | baseline (no batch, obj) | 17075 | 87226 | 1.00× |
  | batch_trace (obj) | 13998 | 85554 | 1.22× |
  | batch_trace+wiring (obj) | 12784 | **16509** | 1.34× |
  | batch_trace (FAST) | 25452 | 85554 | 0.67× ← **S502 loss reproduced** |
  | **batch_trace+wiring (FAST)** | **8873** | **16509** | **1.92×** |

The wiring batch cuts communication **5.3×** (the `K` chain transcripts collapse
to one). The decisive result: globally enabling FAST_BIG is a **loss** with only
the trace batched (the still-tiny per-layer wiring cubes penalise the 24-op
Mersenne mulmod, reproducing the S502 22-vs-16 s loss) and a **1.92× net win**
once **both** kernels are widened — `8.87 s`, well under the 15.5 s baseline. The
S500 prediction ("the right lever is reducing the count of small field-ops, not
the per-multiply cost") is confirmed end-to-end.

**Item-5 headline — sound large-n π(2ⁿ) over `BIG_Q` in the winning config**
(delegate+structured, batch_trace+batch_wiring, FAST_BIG; `claimed π == sieve`):

  | n | x | V=√x | K=π(√x) | wall (s) | π(x) |
  |---|---|---|---|---|---|
  | 16 | 65535 | 255 | 54 | 9.3 | 6542 |
  | 18 | 262143 | 511 | 97 | 22.2 | 23000 |
  | 20 | 1048575 | 1023 | 172 | 66.2 | 82025 |

This is exact π(x) at **x ≈ 10⁶** behind an **Õ(√x) verifier** over a
**sound-characteristic** prime (`BIG_Q=2⁶¹−1` — the wrap-around alias the demo
prime admits above its field is rejected, S498), Õ(x) prover, both wirings
delegated, both big kernels batched. The chain remains prover-bound (the √x-cube
layer sum-checks); the ~3× wall per +2 in `n` tracks the 4× growth in `x`.

**Item-5 reach push to n=22/24 in the FULL succinct config (S513,
`large_x_benchmark.py`).** The table above is the S504 config; the reach driver
runs the **FULL** config — it additionally enables `pcs` (real leaf openings,
S505), `batch_ub` (per-layer Ub-opens → one, S506) and `commit_base` (the two
S₀ base closes via the tensor PCS, S508), i.e. the full non-interactive succinct
certificate (same as `certificate_profile.FULL`) plus `FAST_BIG`. Numbers
(`BIG_Q`, one seed, `claimed π == sieve`):

  | n | x | V=√x | nb | K | π(x) | wall (s, honest) | comm (elems) | comm_outer | per-layer leaf | base opens |
  |---|---|---|---|---|---|---|---|---|---|---|
  | 20 | 1048575 | 1023 | 10 | 172 | 82025 | 71.8 | 73333 | 94% | 0 | 8256 |
  | 22 | 4194303 | 2047 | 11 | 309 | 295947 | 252.8 | 146637 | 95% | 0 | 12352 |
  | 24 | 16777215 | 4095 | 12 | 564 | _(pending background run)_ | | | | | |

The FULL config is ~8% slower than the S504 table at n=20 (71.8 vs 66.2 s — the
pcs/commit_base overhead) but yields the complete succinct certificate. The
profile is the milestone made concrete at the reach: **per-layer verifier
large-table leaf-eval count = 0** (Õ(√x) end-to-end, S506), **comm ≈ 95%
`comm_outer`** (the K sequential reductions, Θ(√x), with every batched discharge
polylog — n=22: `comm_bt/bw/ub` = 80/2373/61), **base opens = the one-time tensor
commitment** (Θ(x^{1/4}), S508). Soundness at the reach is witnessed: the
`delta_pi` liar (claim π+1) is **REJECTED** at n=22 (66.9 s). Details in
`large_x_benchmark_results.md`.

**`commit_code` — the row code for the base commitment (S514).** `run_chain` and
`main` take `commit_code` / `--commit-code` (default `"rs"`, verbatim) selecting
the `pcs_commit` row code: `"rs"` (Reed–Solomon, needs codeword length `N≤q`) or
`"expander"` (Brakedown-style linear-time code, **field-size-free**). The expander
removes the `N≤q` demo constraint that otherwise caps the committable base-table
size; same sub-√x base-open verifier. Chain verdict is UNCHANGED either way —
selftest **case 24** asserts honest `π==sieve`, `delta_pi` + self-consistent liar
rejected, per-layer still leaf-free, base opens still sub-√x
(`vcommit_ops==vleaf_ops_ot`) with `commit_code="expander"` over q AND BIG_Q at
`n∈{8,10,12}`. End-to-end demo: `--n 12 --delegate --structured --field big --pcs
--commit-base --commit-code expander` → `π(4095)=564`, all cheats rejected. Full
expander-code documentation in `pcs_commit_results.md` (S514 section).

## Falsification

Falsified by: the compressed DP disagreeing with a sieve (checked n≤20); an
honest layer **or honest chain** rejected; `claimed π(x) ≠ sieve π(x)`; a
reconstructed output/base claim ≠ the true MLE of `S_{i-1}` (checked both
sides); any cheat — corrupted claim, wiring table, region split, trace
quotient, affine image, small-side gate/division, tampered cross-side carry,
or a wrong claimed π(x)/self-consistent-liar — accepted above the field
soundness bound `~(deg·vars)/q`; or verifier work scaling with `2^n` instead
of `nb`. None observed (selftest covers n∈{8,10,11,12} chains incl. odd-n
padding, two-sided layers, and all cheat classes).

The **field lift** (S499) is additionally falsified by: the chain over `BIG_Q`
(object dtype) giving `claimed π ≠ sieve π` in any of the three modes, or
accepting a chain cheat (it does not — selftest 16); or the lifted prime
**failing to reject** the chain-config wrap-around alias that the too-small
prime accepts (it rejects 5/5 over `BIG_Q`, accepts over `SMALL_Q`); or the
default-`q` transcript differing from the pre-S499 artifacts (it is verbatim —
every prior selftest passes unchanged).

The **verifier op-count curve** (S507, `--bench-verifier-ops`) is falsified by:
config (c) `pcs+batch_trace+batch_ub` having a non-zero per-layer leaf-eval count
(`vleaf_cnt_pl ≠ 0`); its `total vleaf_ops` slope not approaching `0.5` (or the
per-step `Δn=2` growth not `~2×`); configs (a)/(b) not approaching slope `1.0`
(per-step `~4×`); the exact leaf-eval count identities `K·(nb+5)−1` (a), `K·nb`
(b), `0` (c), or the one-time `2·2^nb`, failing (asserted in selftest 22 over `q`
AND `BIG_Q`); the op count differing between `q` and `BIG_Q` (it is field-
independent); `claimed π(x)` differing across the three configs; or the `comm`
slope exceeding `~1.0` (which would re-introduce a Θ(x) non-leaf verifier term —
it is `~0.6`, Õ(√x), in all three).

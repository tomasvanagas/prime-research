# merged_layer_comm.py — results (S510)

## Question

S509 (`certificate_profile.py`) measured the full succinct certificate for
`π(x) = c`: size `α = 0.473` (Θ(√x)), **dominated by `comm_outer`** = the
`K = π(√x)` SEQUENTIAL two-sided layer reductions, after every batchable
obligation (trace/wiring/Ub) was reduced to polylog and the base committed to
`x^{1/4}`. It flagged ONE open quantitative question:

> Is the Õ(√x) certificate size INHERENT to the layering, or just an artifact of
> the unbatched chain? Open problem 2 (CLOSED NEGATIVE) showed merging `j`
> consecutive sieve layers blows up the **prover** (a `j`-layer merge has `2^j`
> fill per row — one entry per subset-product `⌊v/(p_{i₁}···p_{iₛ})⌋`; tree cost
> `m·K²/2 ≈ x^{3/2}`); it never measured the **comm** of a merged-layer verifier.

This cycle measures it as a standalone primitive (one merged layer's comm vs `j`
single layers), not a full-chain rebuild.

## What was built

Work over the SMALL (value-addressed) side, where the sieve wiring is the pure
value→value map `⌊v/p⌋` — the cleanest face of the `⌊v/p⌋` semigroup. The
homogeneous Lucy layer operator is

```
A_p[S](v) = S(v) − [p² ≤ v ≤ V]·S(⌊v/p⌋).
```

(The public Lucy correction `sp_i = S_{i−1}(p_i−1) = i−1` adds only a
verifier-computable additive term per layer — no S-claim — so it does not touch
the `2^j` fill being measured. The homogeneous model isolates exactly that fill;
selftest 1 confirms it equals direct `j`-layer composition.)

Composing `j` layers (`p_1` first … `p_j` last) gives the inclusion-exclusion
expansion with `2^j` subset terms:

```
M[S](v) = Σ_{T ⊆ {1..j}} (−1)^|T| · GateProd_T(v) · S(⌊v/m_T⌋),
   m_T = Π_{k∈T} p_k,
   GateProd_T(v) = Π_{k∈T} g_{p_k}(⌊v / Π_{l∈T, l>k} p_l⌋).
```

Two SOUND reductions of a claim `M̃[S](r) = c` to a claim about the committed
input `S`, both verified against ground truth, both cheat-tested over `q` and
`BIG_Q`:

- **DIRECT** — one degree-3 sum-check over the `v`-cube (`nb` rounds) reduces `c`
  to the `2^j` openings `{ST_T(r')}` (`ST_T[v]=S[⌊v/m_T⌋]`); each is routed back
  to a claim `S(·)` by its own degree-2 wiring sum-check, and the `2^j` S-claims
  are folded (`leaf_open.open_batch`) to one residual closed at the base.
  `comm ~ 2^j · nb`.
- **BATCHED** — stack the `2^j` subset terms on a `j`-bit axis (the
  `batched_trace` pattern): one degree-3 sum-check over the `(j+nb)`-cube reduces
  `c` to a SINGLE opening `STst(R)`, routed to one S-claim. The `2^j` lives
  entirely in verifier/prover WORK (recomputing the public stacked operator MLE,
  `O(2^j·2^nb)`), NOT in comm. `comm ~ j + nb` (polylog).

Three metrics are accounted: `comm` (transcript, split `comm_main`/`comm_fill`),
`vwork` (verifier large-table evals, the S507 metric), and the prover fill-table
size `2^j·2^nb`. The chain extrapolation uses `total = ⌈K/j⌉ · per-merged`. Comm
is **prime-independent and field-independent** (the transcript size depends only
on `nb, j`; selftest 4 asserts identical counts across prime windows and across
`q` vs `BIG_Q`), so the extrapolation is exact, not a representative-group
estimate.

## Headline result (the answer to S509's question)

**The chain TOTAL comm exponent in `x` stays ~0.5 for EVERY merge depth `j` in
BOTH modes — the √x certificate size is LAYERING-INHERENT, on the verification
face exactly as open problem 2 found it on the prover face.**

Fitted exponents `α` (metric ~ `x^α`; `log₂ metric` vs `n`, `n ∈ {10,12,14,16}`,
field `q`; **identical over `BIG_Q`**):

| mode | j | TOTcomm α | TOTvwork α | TOTprover α |
|---|---|---|---|---|
| direct  | 1 | 0.489 | 0.884 | 0.884 |
| direct  | 2 | 0.471 | 0.867 | 0.867 |
| direct  | 4 | 0.470 | 0.867 | 0.867 |
| batched | 1 | 0.481 | 0.884 | 0.884 |
| batched | 2 | 0.457 | 0.867 | 0.867 |
| batched | 4 | 0.446 | 0.867 | 0.867 |

Per-merged growth with `j` at `n=16` (`nb=8`):

| j | 2^j | DIRECT comm_fill | fill/2^j | DIRECT vwork | BATCHED comm_fill | BATCHED vwork |
|---|---|---|---|---|---|---|
| 1 | 2  | 78  | 39.0 | 1792  | 27 | 1536 |
| 2 | 4  | 130 | 32.5 | 3328  | 27 | 2560 |
| 4 | 16 | 442 | 27.6 | 12544 | 27 | 8704 |

- **DIRECT `comm_fill = (2^j + 1)·(3nb+2)`** — clean `2^j` (the inclusion-exclusion
  fill-in). The `fill/2^j` column is flat-ish (the `+1` constant and `comm_main`
  dilute it at tiny `j`; the j-sweep ratios below converge to ×2 per +1 `j`).
- **BATCHED `comm_fill` is FLAT in `j`** (27 everywhere) — the `2^j` left the
  transcript. But **BATCHED `vwork` still tracks `2^j`** (1536→2560→8704, ratio
  5.7× over `j:1→4`) — batching moves the `2^j` out of comm but **not** out of the
  WORK. (DIRECT `vwork` does too: 1792→12544.)

j-sweep at `n=16` (`x=65535`, √x≈2⁸, `K=π(√x)=54`):

| j | groups | DIRECT total comm | BATCHED total comm | BATCHED vwork total | prover fill total | /x |
|---|---|---|---|---|---|---|
| 1 | 54 | 5940  | 3402 | 82944  | 27648  | 0.42 |
| 2 | 27 | 4374  | 1809 | 69120  | 27648  | 0.42 |
| 3 | 18 | 4788  | 1278 | 82944  | 36864  | 0.56 |
| 4 | 14 | 6636  | 1050 | 121856 | 57344  | 0.88 |
| 5 | 11 | 9790  | 869  | 185856 | 90112  | 1.38 |
| 6 | 9  | 15498 | 747  | 299520 | 147456 | 2.25 |

Reading the curve:

- **DIRECT total comm has a shallow minimum at `j=2` then GROWS** — merging
  beyond a constant makes the certificate *bigger* (the `2^j` fill swamps the
  `1/j` layer-count saving). It never drops below the `j=2` value.
- **BATCHED total comm keeps dropping by a constant factor** (3402→747 over
  `j:1→6`) toward the floor `K = 54`: the merge-axis sum-check itself costs
  `log₂(2^j) = j` per merged layer, and `⌈K/j⌉·j = K`. The exponent in `x` is
  still 0.5 — the constant prefactor shrinks but the √x floor does not move.
- **The cost of that comm saving is a `2^j` blow-up in BATCHED verifier WORK and
  in the prover fill-table** (`prover_TOT/x`: 0.42 → 2.25). The prover fill total
  `⌈K/j⌉·2^j·2^nb` trends toward `x^{3/2}` as `2^j → K` — reproducing open
  problem 2's `m·K²/2 ≈ x^{3/2}` on this exact construction.

**Conclusion.** No merge depth, in either batching strategy, drops the comm
exponent below 0.5. Merging trades a bounded comm-prefactor win for a `2^j`
prover/verifier-WORK blow-up. The Õ(√x) certificate is the `K=π(√x)` sequential
`⌊v/p⌋`-semigroup reductions — the SAME wall as open problem 2 (prover-time,
closed negative) and the project's one sequential-sieve wall, now confirmed on
the verification (certificate-size) face.

## Soundness (selftest, over `q = 2³¹−1` AND `BIG_Q = 2⁶¹−1`)

1. the `2^j` subset operator == direct `j`-layer composition (exact ints), `j≤4`;
2. honest merged reduction ACCEPTS, base residual is the true `S̃(r)`, both modes,
   `j≤4`;
3. cheats REJECT in both modes: `c_value` (wrong top claim → sum-check round-1),
   `c_gate` (corrupt a public gate-product table → grounding / round-1 fails),
   `c_route` (self-consistent liar: corrupt one routed `ST_T` and recompute the
   top claim to match → caught by the routing sum-check against the honest `S`);
4. comm prime-independent (different prime windows, same `nb,j`) AND
   field-independent (`q` vs `BIG_Q` identical counts);
5. scaling: total comm exponent in `[0.40, 0.70]` for all `j`, both modes;
   DIRECT `comm_fill` ratio `~2^j`; BATCHED `comm_fill` flat but its `vwork`
   ratio `~2^j`; prover-table exponent ≫ comm exponent.

## What would falsify this

- the subset-sum operator disagreeing with direct `j`-layer composition;
- an honest merged reduction (either mode) rejected, or its residual ≠ true `S̃(r)`;
- any cheating prover accepted above the field bound;
- comm NOT prime-independent at fixed `(nb, j)`;
- DIRECT `comm_fill` NOT growing `~2^j` (e.g. staying linear in `j`);
- BATCHED `comm_fill` growing like `2^j` (it should be flat/polylog);
- **the chain TOTAL comm exponent dropping clearly below 0.5 for some fixed `j` in
  either mode** — THAT would be the genuine breakthrough S509 flagged (a real
  batching of the SEQUENTIAL `⌊v/p⌋` reductions). Verify it against open problem
  2's fill mechanism and re-derive the prover cost before believing it;
- the prover fill-table total NOT trending toward `x^{3/2}` as `j` grows.

## Honest scope

- The homogeneous model drops the public correction scalar (verifier-computable,
  no S-claim) — by design, to isolate the `2^j` S-fill. The full chain's
  per-layer `comm_outer` is `O(nb²)` (line-batched multi-claim + trace/affine
  sum-checks over `Lv` bit tables) vs this model's `O(nb)`; the absolute constant
  differs, but the **j-scaling** (the `2^j` fill) and the **total-comm exponent in
  `x` (0.5)** — the quantities S509's question is about — are the faithful part.
- `vwork` is accounted in the non-delegated model (the verifier recomputes the
  routed/stacked operator MLE directly, `O(2^j·2^nb)`). A delegated wiring would
  cut the per-eval cost but not the `2^j` COUNT of routings — the `2^j` stays.
- Reachable `n ≤ 16` (`nb ≤ 8`, √x ≤ 256), `j ≤ 6` (table build `2^j·2^nb`). The
  exponents are short-sweep fits; the structural formulas
  (`comm_fill_direct = (2^j+1)(3nb+2)`, `comm_batched = 4(j+nb)+const`,
  `total = ⌈K/j⌉·per`, `prover = ⌈K/j⌉·2^j·2^nb`) are exact and match every row.

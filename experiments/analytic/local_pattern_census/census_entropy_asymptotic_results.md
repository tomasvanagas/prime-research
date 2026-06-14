# The census entropy constant via a cluster expansion: a validated estimator, a tight rigorous bracket (S520), and the method's two closure walls (S521)

## Question (open item 4, PROGRAM.md NEXT ACTION after S519)

S519 closed the transfer-matrix route to the *exact* `A_k` (state space
Θ(A_k)) and reduced the count to the active primes `B(k) ~ ρ*(k) ~
k/ln k`:

> `A_k = #{ S ⊆ allowed(k) : for every prime q ∈ B(k), the occupied
> residues {o mod q : o ∈ S} do NOT cover all of Z_q }`.

The remaining question is purely asymptotic — the conjecture
`A_k = e^{(1+o(1))·π(k)}`, i.e. `ln A_k / (k/ln k) → 1`. The S519 NEXT
ACTION asked for an estimate of that constant **without** exact counting,
via a generating-function / cluster-expansion analysis, validated against
the exact `A_k` (k ≤ 64) before extrapolating.

**S520 result.** Built `census_entropy_asymptotic.py`, an exact
reformulation + cluster expansion that (i) reproduces the exact `A_k`
ratio at k ≤ 64, (ii) gives a **rigorous two-sided bracket** on the
entropy ratio that scales past the exact-counting wall, and (iii) shows
the constant's approach to 1 is **driven by the CRT correlations, not the
independent-prime leading term**.

**S521 result (this update).** The S520 NEXT ACTION asked to *tighten the
k=128 bracket toward closure* via (a) a vectorised `joint_cover_prob`, or
*measure whether the bracket width plateaus* (b) — "this directly tests
whether the cluster/block method can EVER prove =1 or whether it provably
cannot." Both are now answered, with a clean **negative** and a
**correction of a prior over-reach**:

1. **(a) done — three engines.** `joint_cover_prob` is now numpy-vectorised
   (≈10× faster, multi-limb uint64 masks for N>64), plus an **exact
   integer** engine `joint_cover_int` (jc·2^N as a big integer, zero
   floating cancellation). All three agree bit-for-bit where comparable
   (selftests [12], [14]).
2. **(b) answered — the bracket width PLATEAUS.** A new `--bracket-sweep`
   measures the best feasible rigorous bracket vs the I-E cost budget `C`.
   The width drops only at **discrete cost thresholds** (when one more
   prime becomes affordable in the bracketing block) and is **flat between
   them**. At k=128 the width is **stuck at 0.0160 for every feasible
   C ≥ 26**; the next improvement needs cost **2³⁹** and the bracket shuts
   only at **2⁷⁵**. The closing cost grows **super-exponentially in k**
   (2²⁶ → 2⁷⁵ → 2²⁷⁹ → 2⁸⁷²). ⇒ **the FKG-block / cluster method cannot
   close the bracket — hence cannot prove the constant = 1 — for any
   k ≥ 128.**
3. **Correction of S520.** A **second, independent wall** was found: at
   k=256 (N=128) the joint miss-probabilities are `~2^{-Θ(N)} ~ 1e-17` but
   the I-E sums terms `~1`, so the **float engine cancels to zero
   significant digits**. The S520 k=256 bracket `[0.9282, 1.126]` was
   therefore **numerical noise, not a rigorous result — it is retracted.**
   The exact-integer engine removes this wall (e.g. exact
   `P_miss({3,5,7})|_{k=256} = 1.48e-19`, where the float gives `0.0`); the
   *corrected* k=256 bracket is below.

## Headline (corrected)

```
   k     ln A_k / (k/ln k)         method                          status
   8        0.6667   (exact)
  16        0.8081   (exact)
  32        0.8861   (exact)
  64        0.9276   (exact)        bracket [0.9236,0.9276] ∋ it
 128     ∈ [0.9243, 0.9403]  RIGOROUS (float OK to ~7 digits at N=64)
                              WIDTH 0.0160, STUCK for all feasible C≥26
 256     ∈ [0.8312, 0.9431]  RIGOROUS (EXACT integer; float underflowed)
                              WIDTH 0.1119; S520's [0.9282,1.126] was NOISE
```

The reliable ratio climbs `0.667 → 0.808 → 0.886 → 0.928`, consistent
with `→1`, but the bracket **does not close** and the constant `= 1` is
**not** proved — and (S521) the cluster/block method **provably cannot**
prove it.

## The exact reformulation (S520, unchanged — the key)

Put the uniform measure on subsets: each allowed offset is included
independently with probability ½. Then

> `A_k = 2^N · P[ admissible ]`,  `N = |allowed(k)| = φ(k)`,  **exactly**.

Let `C_q` = "S covers every class mod q" (*increasing*) and
`\bar C_q` = "S misses a class mod q" (*decreasing*);
`admissible = ∩_{q∈B(k)} \bar C_q`.

- **Marginal:** `p_q = 1 − ∏_j (1 − 2^{−n_{q,j}})`, exact from the class
  sizes (stable via `log1p`/`expm1`).
- **Mean field** = `∏_q p_q`; FKG (decreasing events positively
  associated) ⇒ `P[∩] ≥ ∏_q p_q`, a **rigorous lower bound**.
- **Cluster expansion:** `ln P_B = Σ_{∅≠T⊆B} c_T`,
  `c_T = Σ_{U⊆T}(−1)^{|T|−|U|} ln P_U`; truncating at `|T| ≤ r` is order-r
  (r=1 is the mean field).

## The engine: joint cover by nested inclusion–exclusion (3 implementations)

`P_T` needs the joint cover probs `P[∩_{q∈S} C_q]`, computed by a nested
I-E with the largest prime `L ∈ S` as the product prime; **cost
`2^{(sum S) − max(S)}`** (independent of k):

```
P[cover all q∈S] = Σ_{(U_h ⊆ Z_h)_h}  (−1)^{Σ|U_h|}
                   · 2^{−(#offsets in any head-class U_h)}
                   · ∏_{w∈Z_L} (1 − 2^{−n_w}).
```

- **`_jc_py`** — the readable reference recursion.
- **`_jc_vec`** (S521) — numpy: the largest heads (always including the
  single largest, so a head prime > batch size never falls into a 2^h
  Python loop) are batched into one array of 2^{sum} excluded-offset masks
  via doubling; multi-limb uint64 supports N>64. ≈10× faster (the 2²⁶
  block: 204 s → ~20 s), bit-for-bit equal (selftest [12], rel < 1e-10,
  k=32..256).
- **`joint_cover_int`** (S521) — exact integer `jc·2^N` (each term equals
  `sign·∏_w(2^{n_w}−1)` exactly, since `excl + Σn_w = N`). No floating-point
  cancellation at any N; the I-E numerator equals `count_tm` **exactly**
  (selftest [14]).

Validated bit-for-bit against the S519 transfer matrix `count_tm`
(selftest [2], pairs+triples k=32,64; reldiff 0).

## (b) The bracket WIDTH plateaus — the cost-budget wall (`--bracket-sweep`)

The rigorous bracket at cost budget `C` (= log₂ of the I-E cost cap):

- **Lower** — greedy coarsest partition of `B` into blocks each of cost
  `≤ 2^C`; `P[∩_B] ≥ ∏_blocks P[∩_block]` (FKG). Coarser ⇒ tighter.
- **Upper** — the largest *prefix* block `G` (smallest primes) of cost
  `≤ 2^C`; `A_k ≤ 2^N P[∩_G]` (drop the other primes). Larger ⇒ tighter.

Both advance only when `C` crosses a **prefix-cost threshold** (the cost of
admitting one more prime), so the width is a **step function of C, flat
between thresholds**. Measured at **k=128** (`--bracket-sweep 128`):

```
   C        upper-block    ratio_up   ratio_lo    WIDTH
   8          {3,5,7}        0.9465     0.8267    0.1197
  12          {3,5,7}        0.9465     0.8363    0.1102
  15       {3,5,7,11}        0.9412     0.8888    0.0524
  20       {3,5,7,11}        0.9412     0.8890    0.0523
  24       {3,5,7,11}        0.9412     0.8890    0.0523   (plateau)
  26    {3,5,7,11,13}        0.9403     0.9243    0.0160
  30    {3,5,7,11,13}        0.9403     0.9243    0.0160   (plateau)
  34    {3,5,7,11,13}        0.9403     0.9243    0.0160   (plateau)
```

The upper-block prefix thresholds at k=128 are
`2^{0,3,8,15,26,39,56,75}`. So for **every feasible** budget
`C ∈ [26, 39)` the upper block is pinned at `{3,5,7,11,13}` (ratio
0.9403) and the lower at 0.9243 (the `{17,19,23}` merge needs `2³⁶`,
itself infeasible): **the width is STUCK at 0.0160.** The next tightening
needs `2³⁹` (admit prime 17 to the upper block) and the bracket shuts only
at `2⁷⁵` (the single full-`B` block). The reproduced S520 bracket
`[0.9243, 0.9403]` is therefore the **method's ceiling at k=128**, not a
value that further computation can improve.

**The closing cost is super-exponential in k** (`2^{sum(B)−max(B)}`):

```
   k     |B|   max B   sum B   closing cost   next-prime threshold
  64      5     13       39      2²⁶           (= closing, bracket exact)
 128      8     23       98      2⁷⁵           2³⁹
 256     14     47      326      2²⁷⁹          2³⁹
 512     23     89      961      2⁸⁷²          2³⁹
```

Since `B(k)` grows and `sum(B(k)) = Θ(B(k)²/ln B(k))`, the closing cost is
`2^{Θ(k²/ln²k)}` — **the FKG-block / cluster method cannot close the
bracket (hence cannot prove the constant = 1) at any k ≥ 128.** This is a
clean **closure of the method** — the negative the S520 NEXT ACTION asked
to test for.

The plateau is **engine-independent**: the un-includable large-prime block
`{17,19,23}` (I-E cost `2³⁶`) also blows the transfer-matrix engine — it
aborts past **25.7M states** at k=128 (weak large-prime constraints leave
nearly all subsets alive). Neither engine reaches it; the wall is real, not
an artifact of the I-E.

## The float-precision wall — and the S520 k=256 correction

The joint miss-probabilities shrink as `P[miss within G] ~ 2^{-Θ(N)}`,
while the I-E that computes them sums covering-probs `~1`. The relative
cancellation is `~2^{-Θ(N)}`, so double precision (≈16 digits) holds only
to `N ≈ 60`:

```
   P_miss({3,5,7})     k=64 (N=32)   k=128 (N=64)   k=256 (N=128)
   float I-E           3.69e-4 ✓      3.78e-9 (~7 digits)   0.0   ✗ (underflow)
   exact integer       3.69e-4        3.78e-9               1.48e-19
```

(`reldiff(float,exact)` = 0 at k=64, ~2e-7 at k=128, **1.0 at k=256**.)
So the **k=128 float bracket is reliable** (~7 digits, plenty for a 4-digit
ratio), but **at k=256 the float I-E has zero significant digits**.

**S520 reported a k=256 bracket `[0.9282, 1.126]` "RIGOROUS lower".** That
was computed with the float I-E and is **numerical noise** — `ratio_lo =
1.0856 > 1` is *impossible* for a valid lower bound on a ratio that is
`< 1` and conjectured `→ 1` from below. **It is retracted.** The
exact-integer engine (`--exact-ie`) removes the precision wall; the
corrected k=256 bracket (`--bracket-sweep 256 --exact-ie`, 273 s):

```
   C        upper-block    ratio_up   ratio_lo    WIDTH
   8          {3,5,7}        0.9826     0.6505    0.3321
  12          {3,5,7}        0.9826     0.6638    0.3188
  15       {3,5,7,11}        0.9545     0.7529    0.2015
  20       {3,5,7,11}        0.9545     0.7566    0.1979
  24       {3,5,7,11}        0.9545     0.7566    0.1979   (plateau)
  26    {3,5,7,11,13}        0.9431     0.8312    0.1119   <- best feasible
```

⇒ **corrected k=256 bracket = [0.8312, 0.9431] (rigorous, exact, width
0.1119).** Note the contrast: S520's "rigorous lower 0.9282" is *above*
the true exact lower **0.8312** (and its "upper 1.126" is far above the
true upper **0.9431**) — i.e. S520's k=256 "lower bound" was not even a
valid lower bound, confirming it was float noise. (Like k=128, the width
plateaus: C=20 and C=24 give the same 0.1979; the next change needs 2³⁹.)

Even with exact arithmetic the **cost wall** dominates at k=256: only the
prefix `{3,5,7,11,13}` (5 of 14 active primes) is affordable for the upper
block (`2²⁶`), so the bracket is genuinely **wide** — closing it needs
`2²⁷⁹`. The two walls are independent: exact arithmetic defeats the
precision wall but not the cost wall.

## The structural finding (S520, unchanged and reliable)

The independent-prime (mean-field) ratio **falls** with k (computed from
the *stable* marginals, unaffected by the precision wall):
`0.752 (k=64) → 0.605 (k=128) → 0.394 (k=256)`. So the naive leading term
`∏_q P[miss mod q]` does **not** approach 1 — it declines; the CRT
correlations supply a *growing* fraction of the entropy (≈19% → 35% →
57%). The conjecture's `→1` is **carried by the correlations**, not the
leading term — which is exactly why the bracket is hard to close (the
correlations are `Θ(entropy)`-sized).

## What this delivers / what stays open

- **Delivered (S521):** vectorised (≈10×) + exact-integer joint-cover
  engines; a measured **plateau** of the rigorous bracket width
  (k=128: 0.0160, stuck for all feasible C) ⇒ the cluster/block method
  **provably cannot reach the constant = 1**; the super-exponential
  closing-cost law `2^{Θ(k²/ln²k)}`; the discovery and removal (exact
  integers) of a float-precision wall at N≈100; the **retraction of the
  S520 k=256 bracket** as numerical noise and its **exact-arithmetic
  correction**.
- **Still open (the constant = 1):** unchanged and now *doubly* fenced —
  the cluster/block method cannot prove it (S521 cost wall), so proving
  `= 1` needs an analytic handle on the **full** `Θ(entropy)`-sized
  correlation series as `k → ∞`. This is genuinely new analysis (not new
  code) — the census line is closed for computational cycles.

## What would falsify these results

- The vectorised or exact engine disagreeing with `_jc_py` / `count_tm`
  (selftests [2],[12],[14]; bit-for-bit / exact integer). Not observed.
- The bracket width *not* monotone in the cost budget, or the upper block
  changing *inside* a threshold interval (selftest [13]). Not observed.
- The bracket *not* shutting to the exact value at the closing cost
  (selftest [13], k=32). Holds.
- The float I-E *not* underflowing at k=256, i.e. the precision wall being
  spurious (selftest [14]: raw float `≤ 1e-30` while exact `> 0`). Holds —
  so the S520 k=256 numbers are confirmed noise.
- The FKG bracket failing to contain the exact `A_k` at validated k
  (selftest [10], k=16..64). Holds.

## Files

`census_entropy_asymptotic.py` (one script, `--selftest` (14 groups),
`--k`, `--scan KMIN KMAX`, `--bracket-sweep K [--caps ...] [--exact-ie |
--float-ie]`, `--rmax`, `--cost-cap`, `--no-vec`), this file. Reuses the
S519 `census_transfer_matrix.py` (`allowed_offsets`, `reduce_primes`,
`count_tm`). Cross-refs: `census_transfer_matrix_results.md` (the
`W(B)=ρ*(k)=π(k)` geometry fixing the *exponent*), `local_pattern_census_results.md`.

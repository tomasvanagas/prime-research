# certificate_profile.py — results (S509)

## Question

By S508 the whole-chain **verifier**'s expensive sub-term (large-table `mle_eval`
count) is Õ(x^{1/4}); the verifier side of the milestone is met. The quantity that
is *still* Õ(√x) is the **certificate SIZE** — `stats["comm"]` (the Fiat–Shamir
transcript) has slope ~0.5 in every config (S507/S508 noted this but did not
localize it). This cycle answers: **where, exactly, does the Õ(√x) live, and what
membership statement does it yield for `L_π = {(x,c) : π(x)=c}`?**

## What was built

- `compressed_layer.run_chain` now **partitions `stats["comm"]` by source** (S509,
  purely additive boundary snapshots — transcript and verdict unchanged):
  `comm_outer` (the K = π(√x) sequential per-layer two-sided reductions),
  `comm_bt` (one batched trace zero-test), `comm_bw` (one batched wiring discharge),
  `comm_ub` (one batched Ub-opening discharge), `comm_base` (the two one-time S₀
  tensor-PCS eval proofs). The five sum **exactly** to `stats["comm"]` (asserted on
  every complete certificate).
- `certificate_profile.py` — measures, over n ∈ {8..18}, for the **full succinct
  config** (`delegate+structured+pcs+batch_trace+batch_ub+batch_wiring+commit_base`):
  (i) certificate size in field elements and bits, (ii) verifier large-table op
  count, (iii) prover wall; attributes the comm by source; fits every exponent;
  and contrasts against the S508 "config (d)" (wiring **not** batched).

The full config **strictly extends** the S508 config (d) named in PROGRAM.md by also
turning on `batch_wiring` (built S503/S504 precisely to cut comm). With the wiring
batched, `comm_outer` is *purely* the sequential two-sided layer reductions — the
cleanest form of the result.

## Headline numbers (`python3 certificate_profile.py`, field q = 2³¹−1)

```
  n         x      V  nb    K      pi | cert_elems  cert_KB | verif_ops | prover_s
  8       255     15   4    6      54 |       1337     5.06 |       520 |    0.138
 10      1023     31   5   11     172 |       2777    10.51 |       776 |    0.233
 12      4095     63   6   18     564 |       5259    19.90 |      2064 |    0.460
 14     16383    127   7   31    1900 |      10131    38.34 |      3088 |    0.873
 16     65535    255   8   54    6542 |      18899    71.52 |      4128 |    2.021
 18    262143    511   9   97   23000 |      37983   143.73 |      6176 |    5.269
```

Comm attribution (field elements; `outer` = the K sequential layer reductions):

```
  n     comm |    outer  out% |    bt     bw    ub   base
  8     1337 |      799   60% |    28    344    22    144
 10     2777 |     1911   69% |    36    530    28    272
 12     5259 |     3881   74% |    44    756    34    544
 14    10131 |     8041   79% |    48    949    37   1056
 16    18899 |    16467   87% |    56   1245    43   1088
 18    37983 |    34177   90% |    64   1581    49   2112
```

Fitted exponents (α = d log₂(metric)/dn, last 4 points n = 12,14,16,18; α≈1 ⇒ Θ(x),
α≈0.5 ⇒ Θ(√x), α≈0.25 ⇒ Θ(x^{1/4})):

| metric | α | reading |
|---|---|---|
| **comm (cert size)** | **0.473** | **Θ(√x)** |
| &nbsp;&nbsp;comm_outer | 0.522 | Θ(√x) — **DOMINANT** (60%→90% of comm, →100%) |
| &nbsp;&nbsp;comm_bt | 0.092 | polylog (one batched zero-test) |
| &nbsp;&nbsp;comm_bw | 0.179 | polylog (one batched wiring chain) |
| &nbsp;&nbsp;comm_ub | 0.090 | polylog (one batched Ub opening) |
| &nbsp;&nbsp;comm_base | 0.296 | ≈Θ(x^{1/4}) (tensor-PCS eval proof) |
| verifier ops | 0.258 | Θ(x^{1/4}) — reconfirms S508 |
| prover wall | 0.410 | overhead-bound at reachable n (see scope) |

Effect of batching the wiring (full vs S508 config (d), comm in elems):

```
  n full_comm  full_outer ||   d_comm  d_outer
 16     18899       16467 ||    87944    86757
 18     37983       34177 ||   201015   198790    <- batching wiring: 5.3x smaller comm
```

## Conclusion — the certificate, precisely

**`L_π = {(x,c) : π(x)=c}` has, under a collision-resistant-hash / random-oracle
assumption (Fiat–Shamir of the S491–S508 chain), a NON-INTERACTIVE certificate of
size Õ(√x) field elements (Õ(√x·log²x) bits), whose expensive (large-table)
verifier work is Õ(x^{1/4}).**

The Õ(√x) size is **dominated by — and after batching every batchable obligation,
*is* — the K = π(√x) SEQUENTIAL outer two-sided layer reductions**: `comm_outer`
rises 60%→90% of the total over n=8→18 and →100% asymptotically, with α=0.522≈0.5.
Each reduction is O(nb²)=O(log²x) round scalars, chained K times ⇒ K·log²x =
Θ(√x·polylog). Everything that *could* be batched is polylog (trace/wiring/Ub,
α≤0.18); the base is committed to Θ(x^{1/4}); only the sequential chain remains at
√x.

**Honest sharpening of S508.** The verifier must *read* the Õ(√x) proof — each
transcript scalar is an O(1) sum-check round check — so **total verification is
Θ(√x)**, bounded below by the proof size. The Õ(x^{1/4}) of S508 is the
**large-table-eval sub-term only** (the part that was the milestone's binding term
once the per-layer leaf opens were removed); this cycle makes that explicit. So the
binding Õ(√x) quantity is now the same for size *and* total verification.

**Why polylog is the same wall as open item 2.** Compressing the size below √x to
polylog requires batching the K sequential outer reductions. But each layer's
*output* claim is the next layer's *input* claim — a dependency chain — unlike the
trace/Ub/wiring obligations, which are independent per-layer witnesses sharing one
random trajectory (hence data-parallelizable). Batching the sequential reductions
would be exactly an algebraic compression of the semigroup generated by v↦⌊v/p⌋ —
the layer-batching **closed negative** (open item 2: balanced product-tree fails on
2^j fill-in, cost m·K²/2 ≈ x^{3/2}/ln²x). The certificate-size wall and the
prover-time wall are the *same* sequential-sieve wall, on the two faces of the
problem.

## What would falsify this

- The total `comm` exponent not ≈0.5 (measured 0.473; selftest asserts ∈[0.4,0.65]).
- The per-source attribution not summing to the measured `comm` (asserted exactly on
  every complete certificate; a *rejected* run bails mid-layer and is not a complete
  certificate, so its partial comm is not asserted to close — by design).
- `comm_outer` not the dominant / growing fraction, or its exponent not ≈0.5.
- Any batched-discharge or base exponent ≥ comm_outer's (they must be polylog / x^{1/4}).
- The verifier-op exponent not sub-√x (re-confirms S508; the tight ≈0.25 settles past
  nb≈14 / n≈28 — below that the t-query/reshape constant inflates it, e.g. 0.39 on
  n≤16).
- Claimed π(x) ≠ sieve in any config or field.
- **A claimed batching of the SEQUENTIAL outer reductions that actually drops the
  comm exponent below 0.5** — this would be a genuine breakthrough (it batches the
  ⌊v/p⌋ semigroup), NOT a measurement; verify it against open item 2's mechanism
  before believing it.

## Honest scope

- **Prover wall** measured α=0.41 is *overhead/call-count-bound* (~√x·polylog Python
  round-iteration calls; S500/S501) at reachable n; the asymptotic element-work is
  Õ(x) (established S496/S497), surfacing only past the array-size crossover. The
  certificate **size** is this cycle's headline, not the prover wall.
- The comm **count** (field elements) is field-independent — identical over q=2³¹−1
  and BIG_Q=2⁶¹−1 (selftest §2), differing only by bits/elem (31 vs 61, a constant
  factor; same exponent). Demo runs use q for speed; the sound large-n field is BIG_Q.
- The certificate is **computational** (the base commitment, and Fiat–Shamir, assume
  a CRH/RO). The interactive protocol's per-round soundness is unconditional; the
  Õ(√x) size and the membership statement are stated under the hash assumption,
  consistent with S508.

## Reproduce

```
python3 certificate_profile.py --selftest      # all structural + soundness checks
python3 certificate_profile.py                 # the profile, n=8..18, field q
python3 certificate_profile.py --field big     # over BIG_Q=2^61-1 (61-bit elems)
```

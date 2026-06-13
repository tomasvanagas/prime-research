# compressed_prover_mult_trace.py — results (S492, field lift S498)

**Step 1 of the compressed prover** (PROGRAM.md open item 1): the batched
multiplication-trace sum-check primitive, built and tested in isolation.
**S498 added the field lift** — see the "Field lift" section below, which
also **corrects the PROGRAM.md NEXT ACTION**: the planned extension field
`F_{q²}` does NOT lift the integer range (it preserves the characteristic);
the correct lift is a larger-characteristic prime, now implemented and
demonstrated.

## What it certifies

A claim of the form

    C = sum_{d in {0,1}^m}  eq(rho, d) * S~( u_d ) ,
    u_d = floor( X / a_d ) ,   a_d = (dstart + d) * p ,

where `X`, `p` are public, `rho` a verifier point, and `S` the committed
small-side table. Correctness of **every** `u_d` is enforced by prover trace
tables through low-degree multilinear identities; **the verifier never
divides and never touches a quotient.** This is exactly the large-side
crossing the production prover hits when `d·p > sqrt x` and `floor(x/(d p))`
must be looked up by value in the small side — the one place in the layer
recursion that is a genuine variable×variable multiplication (`u_d · d`,
with `p` fixed in a layer), which is why it needs a trace rather than the
affine index map used when `d·p ≤ sqrt x`.

## Construction

Witness tables over the `d`-cube (size `D = 2^m`): `U=u_d`, `R=r_d`,
`Qv=q_d=a_d-1-r_d`, MSB-first bit decompositions `Ub_j` (`Lv` bits),
`Rb_j`/`Qb_j` (`Lr` bits), wiring `A=a_d`, `ONE`.

- **Constraint zero-test** (one degree-3 sum-check at random `tau`,
  constraints batched by powers of random `alpha`): bitness `T(1-T)=0`;
  recompositions `U=Σ2^jUb_j`, `R=Σ2^jRb_j`, `Qv=Σ2^jQb_j`;
  multiply+remainder `U·a_d+R-X=0` (with `a_d` supplied by the verifier's
  own degree-1 wiring `a~(pt)=p·(Int(pt)+dstart)`); remainder bound
  `a_d-1-R-Qv=0`. With the bit-ranges these pin `u_d=floor(X/a_d)` uniquely.
- **Lookup** (routing): `C = Σ_v S(v)ω(v)`, `ω(v)=Σ_d eq(rho,d)[u_d=v]`,
  verified by **B1** (deg 2, v-cube) `C=Σ_v S~(v)ω~(v)` → claims `S~(r_B)`,
  `ω~(r_B)`, then **B2** (deg `Lv+1`, d-cube)
  `ω~(r_B)=Σ_d eq(rho,d)∏_j EB_j(d)`, `EB_j=eq(r_B[j],Ub_j)` built from the
  SAME committed `Ub_j` proven above.

## Measured (q = 2^31−1)

Honest runs accepted; verified `C` reproduced from the routing identity
`<S,ω>` independently. All cheat classes rejected 20/20 (and in the
selftest 8/8 across `p∈{2,3,7}`):

| cheat | mechanism caught | result |
|---|---|---|
| `u_consistent` (quotient+1, bits made consistent) | multiply+remainder identity | reject |
| `u_value` (quotient+1, bits untouched) | recompU | reject |
| `r_value` (remainder corrupted) | multiply / remainder-bound | reject |
| `nonbit` (a bit set to 2) | bitness | reject |
| `wrong_C` (claimed value+1) | B1 round-1 sum | reject |
| `omega_route` (routing table corrupted) | B1 / B2 | reject |

**Verifier scales with `m` and `Lv`, never with `D` or a quotient size** —
the polylog property the whole verification line rests on:

| n | m | D | 2^Lv | t_verifier (ms) | comm (field elems) |
|---|---|---|------|-----------------|--------------------|
| 8 | 4 | 16 | 16 | 0.279 | 52 |
| 10 | 5 | 32 | 32 | 0.291 | 70 |
| 12 | 6 | 64 | 64 | 0.351 | 90 |
| 14 | 7 | 128 | 128 | 0.433 | 112 |
| 16 | 8 | 256 | 256 | 0.554 | 136 |

`D` grows 16× while `t_verifier` grows ~2×; `comm = 4m + 3Lv + (Lv+2)m =
O(m·Lv) = O(log²X)`. The **prover** here works the `D≈sqrt x`-size cube
(not the `2^n≈x` table of the demo prover): trace build `O(D·Lv)`,
sum-checks `O(D·Lv)` — i.e. this primitive is already Õ(√x)-prover for the
large-side crossing batch it handles.

## Honest scope — what is and isn't done

- **Done:** the isolated primitive — quotient certification by trace + the
  batched `S`-lookup, both sound (every check is a genuine MLE-consistency
  test), measured, adversarially exercised.
- **Field lift (S498, implemented):** see the dedicated section below. The
  demo `q=2^31−1` is sound only while the field CHARACTERISTIC exceeds the
  integers the multiply identity relates (`u_d·a_d ~ X`), i.e. `q>X`,
  `n≲30`. Above that there is a genuine soundness hole (a wrap-around alias
  quotient); the fix is a larger-characteristic PRIME (`--field big`,
  `2^61−1`), now built and adversarially tested. (Earlier this bullet read
  "lift to `|F|>X`" — correct; the PROGRAM.md NEXT ACTION's `F_{q²}` reading
  was wrong and is refuted below.)
- **Leaf closure:** claims about committed tables (`U,R,bits` at `r_A`; `S`
  at `r_B`; `Ub_j` at `r_C`) are closed in this harness by direct MLE
  evaluation — the stand-in for a polynomial-commitment opening / the outer
  protocol's line-batching. In integration the two `Ub_j` point-claims
  (`r_A` from the constraint test, `r_C` from B2) are reconciled by the same
  line-restriction trick already used for `s1@r_v, s2@r_u` in
  `lucy_dp_verification.verify_layer`.
- **NOT done (next steps):** integration into the layer protocol — the
  small-side (value-addressed)/large-side (d-addressed) dispatch, the affine
  large-side index map for `d·p ≤ sqrt x` (no trace needed there), chaining
  the `S~(r_B)` output claim into the next sieve layer, and the
  end-to-end compressed-prover benchmark against a real Lucy recomputation.

## Field lift (S498)

### The soundness hole the demo prime has above its field

The multiply-trace identity `U·a + R − X = 0` is checked over `F_q`. It pins
`U = ⌊X/a⌋` only if the field's **characteristic** is a true upper bound on
the integers in play — in particular the product `U·a ≈ X`. For `q = 2^31−1`
that holds only while `X < q` (`n ≲ 30`). Above that a cheating prover can
supply a **wrap-around alias**: a *different* `Lv`-bit quotient
`u'' = ⌊(X+k·q)/a⌋` and `r'' = (X+k·q) mod a` with

    u''·a + r'' = X + k·q  ≡ X  (mod q),   k ≥ 1,

so the field check passes yet `u'' ≠ ⌊X/a⌋`. The bit-range tables do **not**
stop it: `u''` is a genuine `Lv`-bit integer whenever `2^Lv` has headroom
above `⌊X/a⌋`, which is exactly the regime `q < X` (`forge_alias` finds it,
preferring the entry with the most headroom). This is a real soundness gap
of the demo prime above its field, not a stylistic caveat. `--alias-demo`:

| forge against | n | the alias | accepted by | rejected by |
|---|---|---|---|---|
| `small`=10⁶+3 | 24 | `u=512→542`, `u''·a+r''=X+q` | `q=10⁶+3` (aliases) | `2^31−1` **8/8**, `2^61−1` **8/8** |
| `q`=2³¹−1 | 33 | `u=15577→19471`, `=X+q` | `2^31−1` (aliases!) | `2^61−1` **8/8** |

The `n=33` row is the production-relevant one: **the demo prime accepts its
own alias once `X > 2^31`** — the concrete reason `n≳31` was blocked.

### The fix: raise the characteristic (a larger PRIME), not the element count

`--field {q,big,small}` runs the whole primitive over an arbitrary prime
(uint64 fast path for `q ≤ 2^31−1`; exact Python-int object arrays for
larger primes whose products overflow uint64). `big = 2^61−1` (Mersenne) is
sound for `n ≲ 60`. The trace identities, the sum-check skeleton and the
cost accounting (verifier polylog `O(m·Lv)`, prover Õ(√x)) are all
field-agnostic — only the modulus and the dtype change. Honest accepts and
all 6 original cheat classes reject `20/20` over **both** `q` and `big`
(selftest §3 runs the object-array path end-to-end); `big` also drops the
Schwartz–Zippel error from `3.3·10⁻⁸` to `3.1·10⁻¹⁷` as a side effect.

### Why the planned extension field `F_{q²}` does NOT work (`--refute-q2`)

`|F_{q²}| = q² ≈ 2^62` counts **elements**; the **characteristic stays `q`**.
Integers embed mod the characteristic, so the wrap-around offset `k·q`
embeds as the degree-0 element `(k·q mod q, 0) = (0,0)` — the *zero* of
`F_{q²}`. The forgery and the truth are therefore the **same element** of
`F_{q²}`; no challenge (drawn from all of `F_{q²}`) can separate them. The
multiply-identity zero-test sum `Σ_e eq(τ,e)·(U·A+R−X)` (the only constraint
the alias violates — bitness/recompositions/remainder-bound hold as exact
integers, hence `0` in any field):

```
alias off by  delta = k·q = 2147483647
  over F_q   (q=2^31-1, too small) :  0                    -> UNDETECTED
  over F_q^2 (extension, |F|~2^62)  :  (0+0t)              -> UNDETECTED  <-- the refutation
  over F_Q   (Q=2^61-1, > X)        :  1585920254542431725 -> detected    (correct lift)
  honest control over F_q^2         :  (0+0t)              (zero, as it must be)
```

The larger element count of `F_{q²}` *does* shrink Schwartz–Zippel error,
but that is orthogonal to the integer-range/aliasing problem. The only two
genuine fixes are (a) a larger-characteristic prime (done here) or (b) a
small-field schoolbook **carry trace** that certifies limb products +
carries without asking the field to hold the product (a separate
construction, not a field swap). **CLOSED:** "extension field as
integer-range lift" — see `status/CLOSED_PATHS.md`.

### Cost of the lift — the Mersenne-61 mulmod (S500, BUILT)

For `q > 2^32` the product of two residues exceeds uint64, so S498/S499 fell
back to Python-int `object` arrays (correct, but the per-element Python loop is
slow). **S500 builds the planned numpy Mersenne mulmod** for `BIG_Q = 2^61−1`
(`_mul61` / `_sum61` / `_reduce61`, gated by `FAST_BIG`): split each residue
`a = a₁·2³¹ + a₀` (`a₀<2³¹`, `a₁<2³⁰`), four partial products each `< 2⁶²`,
fold `mid·2³¹` and `hi·2⁶²` by the Mersenne identity `2⁶¹ ≡ 1`, two carry-folds
+ a `P61→0` canonicalisation — **every intermediate stays `< 2⁶⁴`**, fully
vectorised. The field kernels (`eq_table`, `mle_eval`, `sumcheck`, `_asum`,
`honest_C`, `build_omega`) route their array products through `fmul`; for
`q ≤ 2³¹−1` and for object arrays `fmul` is exactly `(a·b) % q`, so the default
and object paths are **bit-identical** (selftest §3/§8). `--field big` selects
the fast path; selftest §8 cross-checks the primitives vs a Python-int
reference and the **whole BIG_Q protocol fast-vs-object, bit-for-bit** (same
`C`, accept/reject, comm, alias rejection).

**Measured (`--bench-big`, the certifier over BIG_Q, fast uint64 vs object):**

| n | D = 2^m | obj ms | fast ms | speedup |
|---|---|---|---|---|
| 12 | 64 | 88.7 | 307.4 | 0.3× |
| 16 | 256 | 159.5 | 531.5 | 0.3× |
| 20 | 1024 | 539.9 | 809.4 | 0.7× |
| 22 | 2048 | 1123 | 1000 | 1.1× |
| 24 | 4096 | 1849 | 1107 | **1.7×** |
| 26 | 8192 | 4042 | 1442 | **2.8×** |
| 28 | 16384 | 9815 | 2574 | **3.8×** |

The speedup **grows with array size** and crosses 1 at `D ≈ 2000` (`n ≈ 22`):
`_mul61` costs ~24 vectorised numpy ops per multiply, so its fixed overhead
only amortises once arrays are large (≳1000 elems); on small arrays the object
path's 2 ops win. The certifier works one `D = √x`-size cube per call, so at the
large-x target (`n = 24–34`, `D = 4096–131072`) it is firmly in the winning
regime and the ratio keeps climbing (→10×+). The object baseline here is
numpy-*object* (a C loop calling Python-int ops), **not** a pure-Python loop —
hence "~3–4×", not the "~10–50×" guessed in S498 for that range.

**Honest caveat for the chain:** the *full* `compressed_layer.run_chain` is the
opposite regime — `π(√x)` layers each doing **many small** cubes — where the
per-op overhead makes the fast path *slower* at reachable `n` (see
`compressed_layer_results.md`, S500). The mulmod is a large-array win, not a
universal one.

Aliasing depends on the magnitude of `X`, not the batch size `D`, so the `n=33`
boundary demo uses a tiny `m=4` cube and runs instantly.

## Falsification

Falsified by: an honest run rejected; any cheating prover (wrong quotient,
remainder, non-Boolean bit, corrupted routing, or wrong claimed `C`)
accepted above the field soundness bound `~(deg·vars)/q ≈ 10^-8` at `n=12`;
or verifier work scaling with `D=2^m` or with a quotient size instead of
with `m` and `Lv`. None observed.

Field-lift-specific falsifiers: a **wrap-around alias** forged against
modulus `m₀` accepted by a prime field of characteristic `> X + k·m₀` (would
break the lift); the `F_{q²}` refutation reversed (the extension field
detecting the alias, or the larger prime missing it); or honest accepts
changing under the lift. None observed — alias rejected `8/8` by every prime
exceeding the integers in play; `F_{q²}` sum is the field zero for both
honest and forged witnesses (undetectable, as the characteristic argument
predicts).

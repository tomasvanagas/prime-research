# pcs_commit.py — tensor multilinear PCS with a sub-√x verifier (S508)

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
- **RS field-size.** RS needs `N ≤ q` distinct points (fine at demo `n`). A
  production system swaps in a field-size-free linear code (Brakedown's expander
  code) behind the SAME `commit/open/verify` interface; the tensor reduction and
  the `O(x^{1/4})` verifier cost are unchanged.
- **Leading-term metric.** `vops`/`vcommit_ops` counts the verifier's dominant
  field-op loops (`O(t(r+k))`); the `O(t·log N)` Merkle hashing is also sub-√x.

## What would falsify this

`verify` disagreeing with `mle_eval` on an honest opening; a forged opening /
wrong claimed value / tampered codeword row accepted above the field bound; the
verifier op-count slope `→ 0.5` instead of `→ 0.25` (i.e. not sub-√x); or
`run_chain(commit_base=True)` changing the chain verdict / `claimed π`.

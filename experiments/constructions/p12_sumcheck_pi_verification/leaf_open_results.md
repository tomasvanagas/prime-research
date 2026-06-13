# leaf_open.py — multilinear-MLE polynomial-commitment OPENING (S505)

## What this is

A real **sum-check opening** for the committed multilinear tables that the
compressed π(x) chain (`compressed_layer.run_chain`) reasons about, built to
replace the `mle_eval(S, point)` **leaf-opening stand-ins** that have stood since
step 1 (the open NEXT ACTION). Two primitives:

- **`open_eval(S, pt, claimed)`** → `(ok, r, residual)`. Proves `S~(pt) == claimed`
  by ONE degree-2 sum-check of the MLE-evaluation identity
  `S~(pt) = Σ_w S[w]·eq~(pt, w)`. It is a **reduction**: it converts the claim at
  `pt` into a polylog transcript PLUS a single **residual** claim `S~(r) == residual`
  at a fresh random point `r`. Verifier work `O(n)` (n = #vars), comm `O(n)`. The
  `O(2^n)` cost moves entirely to the prover's eq-table + round polynomials and to
  the residual (which the caller discharges).

- **`open_batch(S, claims)`** → `(ok, r, residual)`. Folds `k` claims about the
  **same** table into ONE residual via a powers-of-γ random-linear-combination of
  eq tables + one degree-2 sum-check. Drop-in replacement for
  `line_batch_pair`/`batch_on_table`.

- **`close_eval(S, pt, claimed)`** — the trivial `O(2^n)` ground-truth close
  (`mle_eval`), kept for the chain BASE (S₀, a one-time `O(√x)`) and for tests.

## Why a reduction, not a free lunch (honest)

A sum-check opening does **not** by itself make an arbitrary committed-vector
opening cheap: it moves the `O(√x)` work from `pt` to a random `r`. The win is
entirely in WHERE the residual goes:

- **Standalone:** the residual must be closed by `mle_eval` → no asymptotic win
  in isolation (you traded one `√x` eval for a transcript + another `√x` eval).
  This is what the standalone bench measures against — and it already wins as a
  **constant factor** because the opening verifier is `O(n)` while the close is
  `O(2^n)` (see numbers).

- **In the chain (`run_chain(pcs=True)`):** the residual is a claim about the SAME
  committed table at a NEW point — exactly the shape of a carried chain claim — so
  it is **threaded as the next layer's claim** and discharged transitively at the
  chain's base (the two S₀ closes). No commitment, no crypto: soundness rides the
  existing GKR layer sum-checks down to the closed base. This is the "sum-check-
  delegated eq-fold whose own leaf is the next layer's claim" of the NEXT ACTION.

## Measured — standalone primitive

`open_eval` verifier (O(nb)) vs the direct `mle_eval` close it replaces (O(2^nb)),
random tables, `--bench`:

**q = 2³¹−1 (uint64):**

| nb | 2^nb | open verifier µs | mle_eval close µs | speedup | open comm |
|---:|---:|---:|---:|---:|---:|
| 8 | 256 | 8.5 | 55.5 | 6.5× | 24 |
| 10 | 1024 | 9.2 | 77.7 | 8.5× | 30 |
| 12 | 4096 | 10.8 | 126.6 | 11.7× | 36 |
| 14 | 16384 | 14.3 | 277.1 | 19.4× | 42 |
| 16 | 65536 | 14.1 | 1014.0 | 72.1× | 48 |
| 18 | 262144 | 16.7 | 7223.8 | 432.5× | 54 |

**q = 2⁶¹−1 (object dtype, exact Python ints) — `--bench --field big`:**

| nb | 2^nb | open verifier µs | mle_eval close µs | speedup |
|---:|---:|---:|---:|---:|
| 8 | 256 | 9.3 | 198.8 | 21.5× |
| 12 | 4096 | 15.3 | 2301.5 | 150.9× |
| 16 | 65536 | 22.8 | 41410.5 | 1816.8× |
| 18 | 262144 | 52.3 | 188704.3 | 3605.6× |

The opening verifier is **flat in the table size** (the `O(nb)` signature) while
the close is `O(2^nb)`, so the speedup **grows with both nb and the field size** —
exactly the "removes the last √x-/p-linear verifier term" target. Comm is `3·nb`
(n rounds × degree-2+1), tiny and `O(nb)`.

`open_batch` vs sequential `line_batch_pair` folding (`--bench-batch`, q, nb=14):

| k | open verifier µs | open comm | line-fold verifier µs | line-fold comm | speedup |
|---:|---:|---:|---:|---:|---:|
| 2 | 33.1 | 43 | 439.3 | 15 | 13.3× |
| 3 | 38.0 | 43 | 914.0 | 30 | 24.1× |
| 5 | 57.6 | 43 | 1698.7 | 60 | 29.5× |
| 8 | 88.0 | 43 | 2997.1 | 105 | 34.0× |

`open_batch` verifier is `O(k·nb)` (no `2^nb`); the line-restriction fold pays
`(k−1)` closing `mle_evals` = `O(k·2^nb)` on the verifier (and `O(k·nb·2^nb)` on
the prover). Note the honest comm trade: `open_batch` comm is **flat** at `1+3·nb`
in `k`, the line fold's grows with `k`, so they cross around `k≈3`.

## Measured — chain integration (`run_chain(pcs=True)`)

`pcs=True` replaces the per-layer carried-claim leaf openings (the
`line_batch_pair`/`batch_on_table` folds, plus the now-redundant `s2`/`s_B` closes
in `verify_affine_region`/`verify_trace_region`) with `open_batch` whose residuals
thread to the base. Delegated chain, `--bench-pcs`:

| n | V | nb | K | π(x) | tV off ms | tV pcs ms | ratio | comm off | comm pcs |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 15 | 4 | 6 | 54 | 6.87 | 5.55 | 1.24× | 2350 | 2425 |
| 10 | 31 | 5 | 11 | 172 | 19.13 | 15.79 | 1.21× | 6943 | 7108 |
| 12 | 63 | 6 | 18 | 564 | 42.00 | 30.92 | 1.36× | 16144 | 16463 |
| 14 | 127 | 7 | 31 | 1900 | 93.71 | 77.82 | 1.20× | 38115 | 38750 |
| 16 | 255 | 8 | 54 | 6542 | 209.95 | 174.05 | 1.21× | 87226 | 88485 |

End-to-end at n=16 (`--n 16 --delegate --pcs`): π(65535)=6542 verified, t_verifier
174 ms, all 4 chain cheats (delta_pi + self-consistent liar at layers 1/27/54)
rejected 6/6, and the two-sided cheat panel rejected 6/6.

**Honest reading of the ratio.** This is a **stable ~1.2–1.36× constant-factor**
verifier win that does NOT grow with n — because it removes ~5 of the per-layer
direct leaf closes (the 2 carried-claim folds ≈ 3 closes + the 2 redundant
`s2`/`s_B` closes), leaving the dominant per-layer `O(2^nb)` term untouched:
`verify_trace_region`'s **nb Ub-bit-table openings** (line 429). So the per-layer
leaf cost drops from ~`(nb+5)` to ~`nb` direct openings — a real reduction and a
working demonstration of the residual-threading architecture, but **not yet the
asymptotic Õ(√x)**. Comm rises ~1–3% (sum-check folds carry a few more round
scalars than a degree-nb line poly).

## Soundness

- **Wrong claim** (`open_eval`): the round-1 identity `g(0)+g(1)=claimed` fails
  immediately (the honest summand sums to the true value); a prover lying in the
  round polys is caught w.p. ≥ 1−2nb/q (degree 2, nb rounds).
- **Wrong batched claim** (`open_batch`): a wrong `c_i` survives the γ-combine with
  prob ≤ (k−1)/q.
- The eq scalar is **grounded**: the verifier recomputes `eq~(pt, r)` (resp.
  `Σγ^i eq~(pt_i, r)`) in `O(nb)` and checks the prover's folded eq value against
  it, so the eq side of the degree-2 product cannot be faked.
- The residual `scal["S"]` is the prover's IMPLIED `S~(r)` — **not trusted** here;
  the caller discharges it.
- **Chain (pcs):** delta_pi and the self-consistent liar are caught in **phase-A
  round-1** of the relevant layer (the sum over the honest S_{i−1} tables ≠ the
  carried/claimed value) — untouched by the leaf openings — so the verdict is
  preserved; every carried scalar is grounded transitively at the S₀ base.

Tested in `leaf_open.py --selftest` (over q, BIG_Q, SMALL_Q): open agrees with
`mle_eval`; residual equals the true `S~(r)`; wrong claim / corrupted batch claim
rejected; all-zero table edge; agreement on a real compressed-Lucy layer table.
And in `compressed_layer.py --selftest` **case 20**: pcs chain verdict unchanged —
honest == sieve and delta_pi + liar rejected — over q AND BIG_Q, automaton AND
delegated+structured, alone and composed with batch_trace + batch_wiring.

## What this does NOT yet do (the remaining gap)

The single remaining per-layer `O(2^nb)` verifier term after pcs is
`verify_trace_region`'s **nb openings of the Ub bit-tables** (line 429,
`mle_eval(W["tabs"][f"Ub{j}"], r_C)`). These are per-layer **witness** tables, so
their residuals have no carried chain claim to thread to. Removing them requires
folding the Ub openings into the **batched-trace** machinery
(`batched_trace.py`), which already commits all Ub tables across all K layers in
one zero-test — that is the documented next step toward an honestly
**polylog-per-layer**, end-to-end **Õ(√x)** verifier. (Alternatively a real
multilinear commitment — Brakedown/FRI-free — would discharge ALL residuals,
including the S₀ base and the Ub openings, succinctly under a hashing assumption,
trading unconditionality for full succinctness.)

## What would falsify this

`open_eval`/`open_batch` disagreeing with `mle_eval` on an honest prover;
accepting a wrong claimed value or wrong batched claim above the field bound; the
residual ≠ true `S~(r)` for an honest prover; the opening verifier (excluding the
prover's table fold) scaling with `2^nb` instead of `nb`; or `run_chain(pcs=True)`
changing the chain verdict (honest accept / cheat reject) vs `pcs=False`.

## Reproduce

```
python3 leaf_open.py --selftest
python3 leaf_open.py --bench               # open_eval vs mle_eval, vs nb
python3 leaf_open.py --bench --field big   # same over 2^61-1 (object dtype)
python3 leaf_open.py --bench-batch         # open_batch vs line folding, vs k
python3 compressed_layer.py --bench-pcs    # chain verifier: pcs off vs on
python3 compressed_layer.py --n 16 --delegate --pcs   # end-to-end headline
python3 compressed_layer.py --selftest     # case 20: pcs verdict unchanged
```

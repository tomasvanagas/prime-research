# large_x_benchmark — reach results (open item 5)

**Script:** `large_x_benchmark.py` (driver; no new protocol machinery).
**What it runs:** ONE honest verification of `π(2ⁿ−1)` through the FULL
succinct config of `compressed_layer.run_chain`
(`delegate+structured+pcs+batch_trace+batch_ub+batch_wiring+commit_base`)
over the sound-characteristic field `BIG_Q = 2⁶¹−1`, with the uint64
Mersenne fast path (`FAST_BIG`) on; reports wall-clock + the proof's
size/verifier-op profile. Soundness is spot-checked with one `delta_pi`
liar (claim `π(x)+1`) that must be REJECTED.

This is the **reach** deliverable for the "large-x benchmark" line: push
the S504 `n=20` headline (66 s under the batched config) to larger `n`,
recording verified `π(2ⁿ)==sieve` and confirming the asymptotic profile
holds (per-layer verifier leaf-eval count exactly 0 ⇒ Õ(√x) end-to-end;
`comm` dominated by the K=π(√x) sequential outer reductions ⇒ Θ(√x) cert;
large-table opens only the one-time base commitment ⇒ Θ(x^{1/4})).

## Headlines (BIG_Q = 2⁶¹−1, FULL succinct config + FAST_BIG)

| n | x = 2ⁿ−1 | V=⌊√x⌋ | nb | K=π(√x)-ish layers | claimed π(x) | sieve | wall (honest) |
|---|---|---|---|---|---|---|---|
| 20 | 1 048 575 | 1023 | 10 | 172 | 82 025 | 82 025 ✓ | 71.8 s |
| 22 | 4 194 303 | 2047 | 11 | 309 | 295 947 | 295 947 ✓ | **252.8 s** |
| 24 | 16 777 215 | 4095 | 12 | 564 | _(pending background run)_ | | |

All rows: HONEST run **ACCEPTED**, `claimed π(x) == sieve_pi(x)` exactly.
`n=20` reproduces the S504 headline now under the FULL config (S504 used
the same batched bundle without `pcs`/`commit_base`; numbers agree:
π(1048575)=82025).

### Soundness spot-check (n=22)
`delta_pi` liar (prover claims `π(x)+1 = 295948`) is **REJECTED** (66.9 s
— the corrupted top claim fails the chain's base check; the liar run is
faster than honest because the reduction short-circuits on the first
mismatch). Soundness at the reach is therefore witnessed, not assumed;
the full 13-class cheat panel lives in `compressed_layer.py`'s selftest.

## Certificate + verifier profile (n=22, the fully-broken-down row)

```
comm = 146 637 field elems  (~1.15 MB @ 61-bit)
  comm_outer = 139 899  (95%)   K=309 SEQUENTIAL two-sided layer reductions  -> Θ(√x)
  comm_base  =   4 224          tensor-PCS eval proofs for the two S₀ closes  -> Θ(x^{1/4})
  comm_bw    =   2 373          ONE batched wiring discharge (was K chains)   -> polylog
  comm_bt    =      80          ONE batched trace zero-test  (was K)          -> polylog
  comm_ub    =      61          ONE batched Ub-opening discharge (was K·nb)   -> polylog
  (sum = 146 637, exact)

verifier large-table (mle_eval over a √x-size cube) ops:
  per-layer leaf   = 0          <- the Õ(√x) end-to-end property (S506); never grows with K
  one-time leaf    = 12 352
  base-commit opens= 12 352     <- the only x-scaling verifier work, Θ(x^{1/4}) (S508)

reported t_prover  ≈ 8.15 s     (timed regions only; ~3% of wall — wall is field arithmetic
reported t_verifier≈ 1.76 s      outside the timed regions, cf. S501)
```

**Reading.** The headline confirms, at x≈4.2×10⁶ (n=22) and — pending —
x≈1.7×10⁷ (n=24), the three asymptotic claims the protocol stack was
built for:
1. **per-layer verifier leaf-eval count = 0** at every reach n ⇒ the
   per-layer verifier is leaf-free, so the only large-table verifier
   work is one-time (S506); the chain verifier is Õ(√x) end-to-end.
2. **`comm` ≈ 95% `comm_outer`** = the K sequential outer reductions;
   every batchable obligation is squeezed to polylog (`comm_bt/bw/ub` are
   80/2373/61 vs `comm_outer` 139 899). The Θ(√x) cert size is the
   sequential-sieve K layers, exactly as S509/S510 localized it.
3. **large-table opens = the base commitment only** (Θ(x^{1/4}), S508).

## What would falsify this result
- `claimed π(x) != sieve_pi(x)` at any n (the chain mis-counts).
- The honest run REJECTS.
- The `delta_pi` liar is ACCEPTED (soundness broken at the reach n).
- `FAST_BIG` ON changing `claimed`/`comm`/`accepted` vs OFF (selftest
  asserts the fast path is a bit-identical drop-in — checked at n≤12).
- Any per-layer verifier large-table leaf-eval (`vleaf_ops_pl > 0`)
  reappearing (would break the Õ(√x) end-to-end claim).

## Reproduction
```
cd experiments/constructions/p12_sumcheck_pi_verification
python3 large_x_benchmark.py --selftest         # fast correctness/soundness gate
python3 large_x_benchmark.py --n 20             # 72 s   reproduce S504 headline
python3 large_x_benchmark.py --n 22             # ~5 min headline + soundness spot-check
python3 large_x_benchmark.py --n 24 --no-cheat  # ~15 min larger headline, honest only
```
Field default `--field big` (BIG_Q=2⁶¹−1, sound to n≲60); `--no-fast`
forces object dtype (slower at these n; bit-identical verdict).

## Honest scope
- The wall-clock is **prover-bound**, not field- or DP-bound: growth is
  ~3.5×/Δn=2 (n20→n22: 71.8→252.8 s), tracking the π(√x) layers × the
  √x-cube sum-checks. This is the open-item-1 Õ(x) prover working the √x
  state; it is NOT a polylog computation of π(x) (the goal stays blocked).
  The deliverable is the *verification* artifact: a verified exact π(x)
  at x≈10⁷ behind an Õ(√x) verifier / Õ(x^{1/4}) large-table opens.
- `commit_base=True` makes the large-table verifier term computational
  (CRH/RO via sha256, S508); the rest of the chain is unconditional.
- The headline is one seed; the selftest exercises the cheat panel over
  two fields and both fast/object paths at small n.

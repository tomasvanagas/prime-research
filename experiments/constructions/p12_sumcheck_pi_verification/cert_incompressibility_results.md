# cert_incompressibility.py — results

**Question (open item 3, lower-bound face).** S509/S510 proved the chain
certificate's SIZE is Θ(√x) and that this is "layering-inherent" — the K=π(√x)
SEQUENTIAL per-prime reductions (the ⌊v/p⌋-semigroup wall). That is a statement
about the CONSTRUCTION. The complementary question: is a sub-√x certificate
*information-theoretically* impossible, or merely unachieved? A certificate's
SIZE is bounded below by the INFORMATION it conveys, so we separate:

* cert **SIZE** — one transcript block per layer × K=π(√x) ⇒ Θ(√x), regardless
  of information (the construction has K layers);
* cert **INFORMATION** — the joint "hard" bits the transcript must carry that a
  polylog verifier cannot itself predict.

If the information is Θ(√x), no sieve-reconstructing certificate can be sub-√x —
an information barrier. If it is polylog, the √x is purely construction shape and
a sub-√x cert is not ruled out on information grounds (the lower bound would then
have to be COMPUTATIONAL — #P-hardness).

**Proxy.** The chain pins, per layer a=1..K, the large-side survivor count after
sieving by the first a primes — the Legendre partial sieve
`phi(x,a)=#{1<=n<=x : n has no prime factor <= p_a}`. This is the principled
per-layer checkpoint scalar (honest proxy — see scope). We measure the joint
information of `{phi(x,a)}_{a=1..K}` beyond a POLYLOG smooth predictor, at INTEGER
precision (CLOSED_PATHS row 737, S65: relative smoothness does not help — the cert
needs exact values).

Code: `cert_incompressibility.py` (selftest: `--selftest`, 30 checks incl.
Legendre identity, li/R sanity, metric estimators, exponent fitter, synthetic
info-forced-vs-compressible separation, SVD integer-rank, `min_exact_rank`
refactor safety, and density-profile separation of dense vs rank-/bit-concentrated
matrices). Sieve smallest-prime-factor to xmax; sweep x=2^k.

## Headline result (xmax=2^23, k=16..23)

**The √x certificate is INFORMATION-FORCED for any sieve-reconstructing verifier.**
The K=π(√x) checkpoint residuals are integer-INDEPENDENT: even after removing the
best smooth/low-rank model, the residual matrix is FULL RANK once adequately
sampled, so the joint information is Θ(K)=Θ(√x)·polylog.

**Sharpened (S515):** at a wide window (W=192·K) the rank/K **floor is 0.88** with no
downward drift across k=16..22 (rank/K ∈ [0.88, 0.98]) and α_rank/α_K = **0.978 ≈ 1** —
the sampling-robust √x statement is rank = Θ(K) = Θ(√x), independent of the finite-window
log-log discount (α_rank ≈ α_K ≈ 0.40, → 0.5 only as x→∞; §1b). The √x information is also
**DENSE across all K layers** — per-layer hard bits uniform (`bits_uniformity` 0.77, active
fraction 0.98) and prefix integer-rank LINEAR (`rank_half_ratio` 0.52) — so it is carried by
Θ(K) layers, not an o(K) subset (§1c).

### 1. Authoritative measure — integer-reconstruction rank (smooth-model-free)

Build `M[i,a]=phi(X0+i,a)` over a dense window of x (incremental, O(W·K)),
subtract a per-row polylog smooth fit, and find the minimal SVD rank whose
reconstruction rounds to the EXACT integer residuals. SVD subsumes the BEST
low-rank (= best smooth) model, so this needs no analytic Buchstab predictor.
`rank <= min(W,K)`; resolving all K modes needs W ≫ K samples in x.

Window-sensitivity at x=2^19 (K=128) — the residual is full-rank once sampled:

| W/K | 1 | 4 | 16 | 63 | 187 |
|---|---|---|---|---|---|
| rank/K | 0.09 | 0.25 | 0.63 | 0.82 | **0.97** |

Sweep (adaptive window W=64·K, giving rank/K≈0.75–0.97):

| k | x | K=π(√x) | rank | rank/K |
|---|---|---|---|---|
| 16 | 65536 | 55 | 42 | 0.76 |
| 18 | 262144 | 97 | 77 | 0.79 |
| 20 | 1048576 | 172 | 143 | 0.83 |
| 21 | 2097152 | 231 | 223 | 0.97 |
| 22 | 4194304 | 309 | 262 | 0.85 |

**rank exponent α_rank = +0.459 ≈ K exponent α_K = +0.415 ≈ √x (0.5, modulo the
finite-window discount).** rank ≈ K ⇒ the K checkpoint residuals are independent
integer degrees of freedom ⇒ no smaller smooth+low-rank model reconstructs them
exactly ⇒ joint info Θ(K)=Θ(√x)·polylog.

### 1b. Wide-window sharpening (S515) — the floor is the robust √x signal, not α_rank

S511's NEXT ACTION asked to push α_rank from 0.459 toward the ideal 0.5 with a wider
window (W=192·K). Doing so **corrects the premise**: α_rank cannot exceed α_K, and the
**finite-window α_K is itself only ≈0.42** (over k=16..23, `α_K = +0.420`). The reason
is the PNT: K = π(√x) ~ 2√x/ln x, so d log K/d log x = 0.5 − 1/ln x + o(1) — at x=2²², 1/ln x ≈ 0.06,
so α_K reads ~0.42, approaching 0.5 only as x→∞. **The sampling-robust √x statement is
therefore the rank/K FLOOR, not the exponent.** Comparing both windows (xmax=2²³, k=16..22):

| window | α_rank | α_K | α_rank/α_K | rank/K floor..max |
|---|---|---|---|---|
| W=64·K  | +0.459 | +0.415 | 1.106 | 0.75 .. 0.97 |
| W=192·K | +0.396 | +0.405 | **0.978** | **0.88 .. 0.98** |

The wide window **lifts the rank/K floor 0.75→0.88** (no downward drift across the sweep)
and brings **α_rank/α_K from 1.106 → 0.978 ≈ 1**. The 64·K headline α_rank=0.459 was
*inflated* by rank/K climbing across the window (0.75→0.97 from low to high k): rank grew
faster than K only because the low-k points were under-sampled. At W=192·K rank/K is
uniformly ~0.9, so rank tracks K and α_rank ≈ α_K. **Honest reading:** the raw exponent
α_rank = 0.396 at the wide window is BELOW the 64·K 0.459 — this is NOT a weakening but the
removal of the climb artifact. The genuine, sampling-robust result is

> **rank = (0.88–0.98)·K with no decay ⇒ rank = Θ(K) = Θ(π(√x)) = Θ(√x),** and α_rank/α_K → 1
> (both exponents → 0.5 only with the shared, slow PNT 1/ln x discount).

### 1c. Per-layer density — the √x info is DENSE across all K layers, not o(K) of them

The complementary question S511 left: is the joint info spread across the K=π(√x)
checkpoints, or carried by a vanishing subset (e.g. the late high-`a` layers)? Profile of
the residual matrix at x=2²² (W=192·K=59328), `layer_density_profile`:

* **per-layer hard bits ≈ uniform** across the 10 deciles of layer index a/K: 9.8–12.7
  bits/layer, **`bits_uniformity` (min/max) = 0.77**, **active fraction = 0.98** (essentially
  every layer carries ≥1 hard bit);
* **prefix integer-rank grows LINEARLY**: rank(R[:,:⌈fK⌉])/(fK) = 0.84, 0.94, 0.92, 0.90,
  0.88 at f = 0.10, 0.25, 0.50, 0.75, 1.00 — each block of layers adds fresh independent
  directions. **`rank_half_ratio` = rank(first K/2)/rank(K) = 0.52 ≈ 0.5** (dense/linear; a
  value →1 would mean the rank — hence the info — lives in a vanishing prefix);
* residual ENERGY is end-loaded (decile 0 = 0.19, decile 9 = 0.34) — but that is purely the
  scale (φ ~ x at small a, so absolute residuals are larger there). **Bits and rank, the
  scale-free density measures, are uniform**, so the energy loading is not concentration of
  information.

⇒ **DENSE**: the √x is carried by Θ(K) layers, each contributing ~equal hard bits and a
fresh independent rank direction — not a handful of fat layers. This is exactly what makes
the joint info genuinely Θ(K)=Θ(√x).

### 2. Per-x sequence (corroborating, conservative 1D under-sample)

Σ-bits (independent upper bound) α=0.60; **gzip(residual) α=0.388**, AR(3)-residual
α=0.557, empirical-entropy α=0.517 — all track K's α=0.415 (incompressible: the
single K-length sequence does not shrink under gzip/AR/marginal coding). Consistent
with the rank result; weaker because one sequence is a 1-D under-sample of the
(x,a) structure.

### 3. Single-value control — reproduces S36 / SESSION_INSIGHTS(i)

bits(|π(x)−Li(x)|) = 4.3 → 8.0 over k=14..23, slope vs log2 x = **+0.089 ≈ 0**
(value ≈ ½log2 x). A SINGLE π(x) is O(log x) hard bits (|π−R| even smaller,
≤15 over the range). The √x barrier is JOINT across the K checkpoints, **not**
per-value — the S36 reframing ("single value O(log x) bits") is confirmed and the
certificate question is shown to be the genuinely different, joint one.

## Interpretation (precise, with scope)

* **Polylog certificate: RULED OUT** for the sieve-reconstruction class — joint
  info is super-polylog (~x^0.46), so no verifier that reconstructs the sieve
  state can have a polylog certificate. This is the new content: it lifts the S36
  single-value "O(log x) bits, barrier is computational" to the CERTIFICATE
  setting, where the JOINT information is super-polylog.
* **Sub-√x certificate: RULED OUT for this class** — joint info is Θ(√x), so the
  certificate SIZE cannot drop below √x for any sieve-reconstructing verifier.
  This UPGRADES S510: the √x is not only "layering-inherent" (construction) but
  "information-inherent" (the checkpoints carry Θ(√x) bits). Cert-size wall =
  prover-time wall = the sequential-sieve wall, now also from the information side.
* **Single value vs joint:** computing one π(x) is O(log x) bits of information
  yet √x-hard to COMPUTE (S36: a computational, not informational, barrier);
  the CERTIFICATE must convey the joint Θ(√x)-bit state of all K checkpoints,
  which is an INFORMATION barrier. The two are cleanly separated and measured.
* **Relation to S65/row-737:** that row found the φ MATRIX relatively low-rank but
  integer-precision-rank "linear in K" (rank 12 at K=18, ≈35 at K=60). We confirm
  and complete it: with W ≫ K the residual matrix is FULL rank (rank/K→0.97), and
  tie it to the certificate's information floor.

## Honest scope / what this does NOT show

* **Sieve-reconstruction class only.** The bound is on certificates whose verifier
  reconstructs the partial-sieve checkpoints (the chain, and any Meissel/Lucy-style
  recomputation). It does NOT rule out a fundamentally different witness (cf.
  factoring: tiny witness for a hard-to-compute answer). A UNIVERSAL sub-√x lower
  bound still needs the COMPUTATIONAL route (#P-hardness of π(x), open item 3's
  formal half) — exactly the redirection that survives.
* **φ(x,a) is a proxy** for the per-layer transcript. The transcript also carries
  MLE-at-random-point hashes whose information is bounded ABOVE by the state's, so
  H on φ is faithful and, if anything, under-counts redundancy (strengthens the
  "≥√x info" direction).
* **Finite range** (x ≤ 2^23). The finite-window exponents (α_rank=0.40, α_K=0.42 at
  W=192·K) sit below the √x ideal 0.5 because of the PNT 1/ln x log-log discount that
  K=π(√x) itself carries — not a sub-√x rank deficiency (§1b). The sampling-robust √x
  signal is the rank/K floor (0.88, no decay) and α_rank/α_K → 1, both confirmed; the
  exponents → 0.5 only as x→∞.

## Self-correction recorded

An initial run with a FIXED narrow window (W=4096) gave rank ~ x^0.35 and rank/K
drifting down (0.77→0.61), which read as a *sub-√x* info floor. The
window-sensitivity check refuted this: rank is capped by sampling (`rank<=min(W,K)`),
and W=4096 under-samples the K modes at large x. With W∝K the rank is ≈ full. The
sub-√x reading was a measurement artifact; the corrected result is √x-forced.

**S515 correction (the wide-window sharpening).** S511 headlined α_rank = 0.459 (at
W=64·K) and proposed pushing it toward the ideal 0.5 with a wider window. Two things were
slightly off: (i) α_rank cannot exceed α_K, and the finite-window α_K is itself only ≈0.42
(the PNT 1/ln x log-log discount — K=π(√x)~2√x/ln x reads √x as ~x^0.42 at x~2²²), so 0.5
is an asymptote, not a reachable finite-window target; (ii) the 64·K α_rank=0.459 was in
fact *inflated* above α_K=0.415 because rank/K climbed across the window (0.75→0.97) — an
under-sampling artifact at low k, not faster-than-K growth. The wide window (W=192·K) lifts
the rank/K floor to 0.88 (uniform, no drift) and gives the honest α_rank=0.396 ≈ α_K=0.405
(ratio 0.978). So the robust √x statement is the **rank/K floor** (rank=Θ(K)=Θ(√x)), and
the exponent reading ~0.4 is the shared PNT discount — NOT a sub-√x rank deficiency. The
raw exponent went *down* (0.459→0.396), but the result is *strengthened* (floor lifted,
α_rank/α_K → 1, density confirmed).

## Falsifiability

* **The √x floor** is falsified if, at adequate sampling (W ≫ K), the **rank/K floor
  DECAYS with x** (drops toward 0, or α_rank/α_K → 0) — that would mean rank is a sub-√x
  power of K and a sub-√x cert is not ruled out. REFUTED here: at W=192·K the floor is
  0.88 with no downward drift and α_rank/α_K = 0.978 ≈ 1. **NB (S515 correction):** the
  raw α_rank crossing below 0.45 is NOT a valid falsifier — α_K itself is ≈0.42 over the
  reachable window (PNT 1/ln x discount), so α_rank ≈ 0.40 at the wide window is the
  *expected* tracking value, not a weakened floor. The correct test is on rank/K decay
  and the α_rank/α_K ratio.
* **Density** is falsified if the √x were carried by an o(K) subset of layers: a
  `rank_half_ratio` → 1 (rank lives in the prefix) or `bits_uniformity`/`active_frac` → 0
  (a few fat layers carry all bits). NOT observed: `rank_half_ratio` = 0.52, `bits_uniformity`
  = 0.77, `active_frac` = 0.98 (uniform across deciles). The selftest checks both metrics
  separate a dense matrix from rank-concentrated (duplicate columns) and bit-concentrated
  (sparse) controls.
* If the per-x gzip / AR exponents were ≈ 0 (residual compressible to polylog),
  the √x would be construction shape only — NOT observed (gzip α=0.39, AR α=0.56).
* If a different (non-sieve-reconstructing) witness for π(x)=c with sub-√x size
  existed, this measurement would be silent on it — that is the explicit scope
  boundary and the reason the #P-hardness route is the remaining universal lever.
* Reproduce: `python3 cert_incompressibility.py --xmax 8388608` (full sweep +
  window-sensitivity + 64·K & 192·K rank sweeps + per-layer density profile, ~3 min);
  `--wfactors 192` for the wide window alone; `--selftest` for the 30 unit checks.

# Critique — Session 49 Round 2 (Proposals A–C, TT/LASSO/learned-residue)

Critic-mode evaluation of the three NEW Session 49 proposals (proposals_session.md):
- **A** — Reordered tensor-train compression of `delta(n)` across 7 orderings;
- **B** — LASSO compressed-sensing of `delta(n)` in the zero-mode basis;
- **C** — Learned (ridge) residues at fixed zeta zeros.

The proposer self-reported all three CLOSED with empirical numbers. This
critique:

1. cross-checks each against `status/CLOSED_PATHS.md` (705 lines / 533+ entries),
   citing the exact prior closures (most are pre-Session-49);
2. flags **Proposal C's basis misspecification** (cos(γ log n)/√n uses prime
   *index* n, but the explicit-formula basis is in prime *value* p(n));
3. runs the corrected basis and confirms the negative verdict survives;
4. preserves the bond-profile saturation pattern and LASSO `nnz=0` finding
   so future sessions can reference them;
5. adds three CLOSED_PATHS entries that pin the new measurements.

---

## A — Reordered TT compression of delta(n)

**Verdict: DUPLICATE — the multi-ordering test is a refinement, not a new
mechanism. All seven orderings produce identical volume-law saturation,
confirming line 171's S10 closure is ordering-independent.**

**CLOSED_PATHS adjacency:**

| Line | Entry | Match |
|------|-------|-------|
| 171  | MPS / tensor network — FAIL/I; bond dim ~ N^{0.49}, volume-law (S10) | EXACT mechanism (natural ordering only) |
| 185  | Tensor sieve (MPS of divisibility DFAs) — FAIL/I (S20) | Adjacent: DFA-product variant |
| 516  | Tensor network / MPS sieve — FAIL/I; bond = primorial (S20) | Adjacent: sieve-not-delta variant |
| 517  | MPS bond dim of chi_P at base-W half cut — FAIL/E; theorem gives rank = min(W^j, φ(W)·W^{d-j-1}+1) (S41) | Theorem-level closure for any base |
| 599  | Tensor network MPS Legendre sieve contraction — FAIL/I; bond ~ 2^{0.33-0.43·a} (S26) | Adjacent: sieve-state target |

**The new measurement.** At N=2^13=8192, eps_rel=1e-3 SVD truncation, **all
seven orderings** (identity, bit-reversal, Gray, Morton, 2-adic, sort_by_R^{-1},
random) plus Gaussian noise produce **identical** bond profiles
`[2, 4, 8, 16, 32, 64, 64, 32, 16, 8, 4, 2]`. The smooth-cosine positive
control achieves bond-dim 2 throughout, confirming the SVD pipeline is
sensitive to structure when present.

**Why the new test does not break the closure.** The bond profile
`[2, 4, 8, 16, 32, 64, 64, 32, …, 2]` is the *mathematical maximum* at each
cut for a binary L=13 TT (cut k bonded by min(2^k, 2^{L-k})). Saturating
this at threshold 1e-3 is necessary but not sufficient for volume-law; the
gaussian-noise baseline saturates the same ceiling, indicating delta is
indistinguishable from GUE noise at this resolution. Tighter thresholds
(1e-6, 1e-9) would still fail to expose structure: line 24 (S25) shows
DFT spectral flatness 0.93–0.999 on zeros, so any reorthogonalization of
delta inherits a flat singular-value spectrum.

**Independent re-derivation.** Theorem in line 517 establishes the chi_P
bond-dim at every cut for any base-W half-cut; the new test is the same
formula instantiated for W=2 with index permutations, all of which give
the W=2 ceiling 2^min(j, L-j). Bit-reversal/Gray/Morton are 2-adic
isomorphisms of the binary cube and cannot change the half-cut min;
sort_by_R^{-1} reorders R^{-1}(n) into a near-monotone curve but the
*residual* delta = p(n) - R^{-1}(n) is GUE-random regardless of how the
indices are sorted. The 7-ordering test cannot succeed because the
residual after smooth subtraction is invariant under any
order-preserving sort applied to the sum, and bit-permutation orderings
permute *labels* without changing the SVD spectrum of the multilinear
contraction.

**Useful refinement to keep.** The seven-ordering identical-saturation
result is the most explicit statement of "volume-law is ordering-independent
under binary TT" the project has on record. **Bond profile
`[2,4,8,16,32,64,64,32,16,8,4,2]` is the binary-TT volume-law fingerprint**
across N=8192 and all seven tested orderings.

**Failure mode:** I (information-theoretic).

**Action:** add CLOSED_PATHS entry sharpening line 171 with the seven-ordering
identical-saturation result.

---

## B — LASSO compressed-sensing on zero-mode basis

**Verdict: DUPLICATE — direct hit on TWO existing CLOSED_PATHS entries.**

**CLOSED_PATHS adjacency (in age order):**

| Line | Entry | Match |
|------|-------|-------|
| 350  | Compressed sensing on K zeros — PARTIAL/I; "lucky cancellations small x, doesn't scale" (S7) | EXACT mechanism, original closure |
| **699** | **Compressed-sensing recovery of psi(x) via random K-subset of zeta zeros** — FAIL/E; first-K is near-optimal among K-subsets, signal DENSE in zero basis with smoothly-decaying amplitudes (S54) | **EXACT mechanism, last-session closure** |
| 623  | Compressed sensing on prime indicator — FAIL/C; "99% energy in 1.25% Fourier components but those ARE zeta zero frequencies. Circular." (S29) | Adjacent: prime-indicator target |
| 472  | Goldreich-Levin / list decoding for pi(x) — FAIL/I; commrank 2^{N/2} blocks any sparse Fourier recovery (S18) | Structural: blocks any sparse-basis route |
| 240  | Wavelet/Fourier compression of correction C(x) — FAIL/I; 99.9% energy needs N^{0.75} coefs (S40) | Adjacent: same target signal |
| 653  | Contour integral evaluation of zero sum — FAIL/E; "99% energy needs 133/500 components. No sparse representation in any basis." (S32) | Adjacent: same basis, different fitness |

**The new measurement.** At N=4096, K_max=1024, alpha sweep
{1e-6, …, 1e-1}, the LASSO best-test-RMSE solution at every K and every
alpha is **nnz=0** (predict delta=0). OLS overfits catastrophically at
K≥256 (test RMSE blows up to 1e10). The naive `round(R^{-1})` baseline
matches the all-zero LASSO at 1.1% prime recovery.

**Independent verification of the proposer's negative.**
Line 699 last session ran a closely-related test: random K=50 subsets of
2000 zeros for psi(x) recovery, finding random-subsets STRICTLY WORSE than
first-K (rms 404 vs 97, 4× worse). Conclusion: signal is dense in zero
basis with 1/|rho_k|-decaying amplitudes, **not K-sparse**. LASSO on the
same basis with no amplitude prior simply discovers the same fact in a
weaker form: every coefficient must scale ~1/|rho_k| to match the
analytic explicit formula, so L1-promoting sparsity is unable to prefer
any small subset.

The sklearn LASSO at alpha=1e-1 is sufficiently aggressive to zero
everything. At alpha=1e-6, the model approaches OLS and overfits as soon
as K > N_train/2 = 1024. There is no mid-alpha sweet spot because there
is no actual sparsity to discover.

**Failure mode:** I (information-theoretic) reinforced by E (signal has no
sparse basis representation by line 699 evidence).

**Action:** add CLOSED_PATHS entry refining line 350 with the LASSO
all-zero outcome and line 699 cross-reference.

---

## C — Learned (ridge) residues at fixed zeta zeros

**Verdict: DUPLICATE + FLAWED BASIS — the proposer's basis
`cos(γ log n)/√n` uses the prime *index* n, but the explicit-formula's
natural variable is the prime *value* p(n). The qualitative conclusion
survives the correction (verified empirically below), but the experiment
as written does not test what it claims.**

### The basis-misspecification flaw

The explicit formula gives
`f(x) = π(x) − R(x) ≈ −2√x/log x · Σ_k cos(γ_k log x − φ_k)/|ρ_k|`.

Inverting around the prime sequence: for `δ(n) = p(n) − R^{-1}(n)`,
linearizing R^{-1} at n gives
`δ(n) ≈ f(p_n) · log(p_n) ≈ −2√p_n · Σ_k cos(γ_k log p_n − φ_k)/|ρ_k|`.

The natural ridge basis is therefore **`{√p_n · cos(γ_k log p_n),
√p_n · sin(γ_k log p_n)}`** — uses prime *value* p_n (not index n) and
amplitude `√p_n` (not `1/√n`). The proposer used `cos(γ log n)/√n`, which
is the *wrong* basis: amplitude scales like `1/√n` ≈ `1/√(N/log N)`
instead of the correct `√p_n` ≈ `√(N log N)` (a factor of N off).

### Critic experiment with corrected basis

`experiments/proposals/critique49_basis_misspec.py` runs ridge regression
under both bases at K ∈ {32, 128, 512, 2000} on the same train/test split.

| K | orig acc | correct acc | orig test RMSE | correct test RMSE |
|---|---|---|---|---|
| 32 | 0.015 | 0.014 | 37.40 | **25.32** |
| 128 | 0.011 | 0.017 | 37.39 | **23.66** |
| 512 | 0.012 | 0.000 | 38.12 | 1670.74 |
| 2000 | 0.011 | 0.000 | 38.27 | 795.42 |

Naive baseline: 0.011.

**Findings:**
- Correct basis gives strictly lower test RMSE at K=32, 128 (25.3 vs 37.4;
  23.7 vs 37.4) — the correct basis IS more informative.
- Correct basis catastrophically overfits at K ≥ 512 — same overfitting
  pattern as OLS, regularization cannot save it.
- **Even the most-informative correct basis at K=128 gives 1.7% recovery,
  ≈naive 1.1%** — RMSE 23.7 means predictions are off by ±24 on average,
  far above the 0.5 rounding threshold.

**The basis misspecification did NOT hide a positive result.** Proposal C
closes under both the original and the corrected basis.

### CLOSED_PATHS adjacency for the underlying mechanism

| Line | Entry | Match |
|------|-------|-------|
| 145  | Deep ML (ridge/AR) — FAIL/I; 5.4% test exact (S8) | Direct mechanism: ML on delta |
| 561  | ML (Ridge/kNN) prediction of delta(n) — FAIL/I; 0% exact on 1000 test cases (S21) | EXACT: ridge on delta |
| 622  | Polynomial empirical correction to R(x) — FAIL/I; "Polynomial in 1/log(x) reduces error 60-80% but cross-validation confirms overfitting" (S29) | Adjacent: learned coefficients on a fixed basis |
| 664  | Spectral truncation / adaptive zero selection — FAIL/E; "1000 zeros insufficient for n≥500" (S33) | Adjacent: same basis, no learned coefs |
| **691** | **Extended PSLQ basis for f(x) — 10 zeta-zero modes** — FAIL/I; cross-check residuals 4e2-6e3 (S48) | EXACT basis structure |
| **692** | **Hermite/Gaussian/Riesz mollification of explicit formula** — FAIL/E; mollification at most 1.21x better (S48) | Direct: re-weighted-coefficient variant |
| **702** | **PSLQ + structured zeta-zero/log basis on delta(n) (8th PSLQ-on-delta)** — FAIL/I+E (S49 round 1) | Direct: PSLQ-coefficient variant of same atoms |

**Why no learned-coefficient choice can break the bound (re-derivation).**
The truncation tail beyond K zeros has L^2 norm
`≥ Σ_{k>K} 1/(γ_k^2) ~ 1/(K log K)`. Multiplying by the `√p_n` amplitude
gives an irreducible test-RMSE floor of `Ω(√(p_n)/(√K · √log K))`. At
N=4096, p_N≈40000, K=2000, this floor is `~ √40000/√(2000·7.6) ~ 1.6`, so
ridge plateaus near naïve no matter what coefficients are learned within
the kept-K subspace.

The "learned coefficients can absorb tail energy" hypothesis is
**provably false in L^2**: the tail energy is orthogonal to the kept-K
subspace, so the L^2-optimal kept-K coefficients are exactly the analytic
1/ρ residues, and any data-fit deviates further from the truth on the
test window. This is why the proposer found **at most 0.4pp** gap
between analytic and ridge: that gap is finite-sample noise.

**Failure mode:** I (information-theoretic) — same as line 691, 692.

**Action:** add CLOSED_PATHS entry refining lines 691 and 692 with the
ridge-vs-analytic null gap and the basis-misspecification correction.

---

## Summary

| Proposal | Verdict   | Failure Mode | Strict prior closure | New CLOSED_PATHS entry? |
|----------|-----------|--------------|----------------------|--------------------------|
| **A — Reordered TT** | DUPLICATE | I | line 171 (S10) | YES (refines: 7 orderings identically saturate volume-law ceiling at N=8192) |
| **B — LASSO on zero modes** | DUPLICATE | I | lines 350 (S7), **699 (S54)** | YES (refines: LASSO best-test-RMSE = nnz=0 at every K, alpha) |
| **C — Learned residues** | DUPLICATE + FLAWED BASIS | I | lines 561 (S21), 691/692 (S48), 702 (S49 R1) | YES (refines: corrected `√p_n · cos(γ log p_n)` basis still plateaus at naïve) |

**No proposal survived critique with a polylog or sub-√x route to π(x).**
All three close to entries that pre-date this session — Proposal B's mechanism
was last closed *last session* (line 699, S54), suggesting the proposer-mode
prompt failed to surface the most recent CLOSED_PATHS entries.

**Process note.** Proposal A's seven-ordering test is the strongest empirical
"ordering-independent volume-law" demonstration on record — preserve as a
refinement. Proposal B's LASSO discovers L1-sparsity is impossible by getting
nnz=0 — duplicate of S54 line 699 but worth recording. Proposal C contains
a basis-misspecification flaw (uses prime index n instead of prime value
p_n); critic re-ran the corrected basis and confirms the closure stands.

**Open status unchanged.** The only genuinely open avenue per
`status/OPEN_PROBLEMS.md` is **circuit complexity of π(x)** (Berry–Keating /
TC^0 lower-bound monitoring). All three proposals attack the analytic /
information-loss side; none addresses circuit complexity.

**Process recommendation.** The proposer wrote `cos(γ log n)/√n` from
memory, not from the explicit-formula derivation. A pre-flight check that
"the basis used matches the basis derived from the explicit formula" would
have caught this. The result happens to be correct (closure verified by
critic experiment), but on a future iteration the misspecification could
hide a positive signal in a different proposal.

## Files written this session

- `archive/ephemeral/critique_latest.md` (this file).
- `experiments/proposals/critique49_basis_misspec.py` + `_results.md`
  (corrected-basis ridge test for Proposal C).
- Three entries appended to `status/CLOSED_PATHS.md`.

## Files NOT written

- No new `novel/` files — every proposal is a strict duplicate.
- No new `experiments/proposals/session49_*.py` — proposer's three
  experiments adequately closed each path.

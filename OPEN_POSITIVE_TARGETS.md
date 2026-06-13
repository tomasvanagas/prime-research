# Open Positive Targets — π-Related Problem Variants Where Algorithmic Wins Are Plausible

After 200+ sessions and five closed commit threads, the headline
polylog π(x) goal is provably blocked along every path the project has
explored. **But the Correlation Dichotomy (S224, Thread 5) showed that
adjacent problems CAN admit partial-positive results** — a 33× speedup
for batched correlated π queries, conditional batch-polylog for narrow
windows.

This file catalogues additional adjacent problems where similar
partial-positives are plausible. Each entry is shaped as a candidate
Correlation-Dichotomy-style attack: a problem distinct from single-x
polylog π(x), with a falsifiable empirical-first attack plan, and a
measurable A-grade outcome.

**These are the framework's current highest-EV positive-direction
targets.** When `commit` mode is invoked without an active thread, or
when `reframe` mode fires, agents pick from this list.

---

## §1. Batched & Windowed Variants

### P1 — π in arithmetic progressions, batched on modulus q [CLOSED-B at S231]
**Problem:** for fixed x and a fixed family of moduli `Q = {q_1, ..., q_M}`
with q_i ≤ x^c, compute `π(x; q_i, a_i)` for all i. Single-(q, a)
queries are well-studied; *batched* on Q is not.

**Status: CLOSED-B at S231 (Thread 6 final wrap).** P1 closes
negatively across all four amortisation axes:

- **AXIS 1 (a-direction, S227):** trivial 8× at M=256 — zeros
  independent of a, predicted exactly by E1.5.
- **AXIS 2 (χ-direction fixed q, S228+S229):** bounded constant
  1.79-2.04× via cyclic DFT primitive over (Z/qZ)*. BLAS-vs-Bluestein
  FFT FLOP rate sets the ceiling.
- **AXIS 3 (composite q, S230):** bounded constant 1.75× at q=1001
  via CRT-based multi-axis FFT primitive over direct-product groups.
  Same crossover as cyclic FFT; no asymptotic gain.
- **AXIS 4 (Q-batched cross-conductor, S231):** bounded constant
  decreasing toward 1× as M grows (1.25× → 1.12× → 1.05× at
  M=2,3,4). Stages B-G all per-conductor; saving from sharing
  Stages A+E is O(1) per query. saving/total → 0 as M·√q_avg grows.
  OPPOSITE shape from S224 Correlation Dichotomy.

**Why P1 fails where P5 (Thread 5) succeeded:** Thread 5 operates on
SHARED-L geometry (same ζ-zero database across queries). P1 / Thread 6
fails on DISTINCT-L geometry (distinct L(s, χ_q) per query, no shared
zero database across distinct characters or conductors).

**Original A-grade target (NOT achieved):** per-(q, a) amortised polylog
over M = poly(log x) characters. The five slots produced bounded-
constant lifts only (1.05-2.04×).

**5-session commit thread spent. See `archive/sessions/session{227,228,229,230,231}_commit_p1_*.md`
for slot-by-slot synthesis. Thread 6 marked DONE_PARTIAL_NEGATIVE in
.commit_state.**

### P2 — Prime gap function `π_h(x) = #{p ≤ x : p+h prime}` for many h [Thread 8 CLOSED-CONDITIONAL at S421]
**Status: CLOSED-CONDITIONAL — Thread 8 done at S421, all 4 slots
complete. Aggregate result: polylog-time HL evaluator HL_∞(h, x) =
S_h · li_2(x) for h-batches H with max h ≤ (log x)^β, |H| ≤ poly(log x),
having cross-h L²-typical error ε_typ(x, H) ≤ F_H · √(S̄_H · x) / log x,
**conditional on the cross-h Hardy–Littlewood Random-Residual Hypothesis
(HLRH)**. Empirically verified across two decades (S419 / S420), with
slot-3 named-exponent variance decomposition σ²_HL_Q = σ²_∞ + ⟨ε_Q²⟩ ·
li_2² rigorised under HLRH, and the knee Q* ≍ √x/log x derived from
quadrature equality. See `experiments/analytic/batched_pi_h/
slot4_theorem.md` for the conditional theorem statement and proof
outline.**

**Problem:** for fixed x and many h ∈ {2, 4, 6, ..., H}, compute
`π_h(x)` simultaneously.

**Why partial-positive plausible:** the singular series factor
`S_h = 2 C_2 Π_{p|h, p>2} (p-1)/(p-2)` for distinct h shares small-prime
factor structure. A sieve that produces all `π_h(x)` for h ≤ H may
share the underlying sieve table efficiently.

**Slot 1 (S418) result — structural dichotomy.** Empirical baseline at
x ∈ {10⁵, 10⁶, 10⁷}, H ∈ {20, 50, 100, 200, 400} (plus H = ⌊log²x⌋
per-anchor) shows a **clean two-regime split**:
- **EXACT regime** (batched-sieve): per-h amortised cost
  T_batched/M = 0.6 ms → 4.6 ms → 38.9 ms across x = 10⁵ → 10⁶ → 10⁷
  with M = 66, 95, 100. Ratio matches π(x) growth (8.18×, 8.47×).
  H-sweep at x=10⁶ shows the per-h amortised cost is **M-independent
  for M ≥ ~50** (8.7→6.0→5.0→4.5→4.3 ms at M = 10, 25, 50, 100, 200).
  M1/M2 batched-vs-per-h speedup is bounded constant (6.7× and 8.85×),
  not growing with M. **Per-h amortised cost is Θ(x/log x), not
  polylog.** This is the P1 (Thread 6) negative shape transposed
  to the h-axis.
- **APPROXIMATE regime** (HL = S_h · li_2(x)): per-h cost
  0.5–1.2 µs, **flat in x and M** — the polylog shape. Mean and max
  approximation errors |π_h − HL_h|/√x ≤ 0.10 and ≤ 0.25 across all
  261 (x, H, h) cells measured. Matches Hardy-Littlewood's O(√x)
  conjectural bound.

**Structural conclusion (S418):** P2 is **Thread-7 (P3) shape, not
Thread-5 shape** — the positive-direction win is the polylog
HL-approximation, with √x error trading away log²x bits of accuracy.
The cross-h amortisation in the EXACT regime gives only a
constant-factor sieve-sharing lift (matching P1).

**Original A-grade target (NOT achievable for EXACT):** per-h amortised
polylog. Slot 1 establishes a Θ(x/log x) lower bound for the natural
batched-sieve algorithm; non-sieve / quantum / GPY-weight approaches
are not yet ruled out but are outside slot-1 scope.

**Pivoted A-grade target (ACHIEVABLE for APPROX, in progress):** a
precise named-exponent error bound for the HL approximation,
analogue of Thread 7 Corollary B but for π_h(x).

**Slot 2 (S419) result — HL residual distribution shape.** Multi-anchor
multi-h ensemble at x ∈ {10⁶, 10⁷, 10⁸} × h ∈ {2, 6, 14, 22, 30, 42,
90, 198}, N=30 log-uniform x-samples per anchor in a quarter-decade
window. Three statistical aggregates:
- **Cross-x at fixed (anchor, h) FAILS half-Gaussian KS:** median p_eff
  = 0.033, only 8/24 cells with p > 0.1. Diagnosis: residuals are
  random-walk-like correlated within a quarter-decade window (range
  metric (max|r|-min|r|)/σ_eff ≈ 1.5 mean across cells).
- **Cross-h at fixed (anchor, x_j) PASSES half-Gaussian KS:** median
  p_eff = 0.7-0.8, **all 90/90 cells with p > 0.1** across three anchors.
  This is the natural ensemble for HL residual analysis.
- **σ_eff / σ_pois ratio in [0.36, 0.70]** across (anchor, h) cells.
  Empirical RMS is consistently SMALLER than HL random-residual
  Poisson prediction √(S_h · li_2(x)) — analogue of Thread 7's GUE
  pair-correlation suppression. NOT stable across decades (0.50 →
  0.70 → 0.41 anchored at 10⁶/10⁷/10⁸; the 10⁷ outlier driven by a
  large h=6 cousin-prime residual episode).

**Structural conclusion (S419):** the right ensemble for HL residual
on the h-axis is **cross-h at fixed x**, not cross-x at fixed h. Half-
Gaussian shape holds on cross-h; the σ-suppression factor analogue of
F_GUE = 0.755 is in [0.36, 0.70] but lacks Thread 7's flat decade-
stability. Thread-7-shape A-grade target (named-exponent error bound)
remains achievable; slot 3 should target the Q-truncation cost-vs-error
curve as the path there.

**Slot 3 (S420) result — Q-truncation tradeoff curve.** 26-h ensemble
spanning max-prime-factor [3, 2003] × Q ∈ {30, 50, 100, 200, 500, 1000,
2000, 5000, ∞} × x ∈ {10⁷, 10⁸} (5 x-samples per anchor). Output:
- **Named-exponent variance decomposition** σ²_HL_Q(x) = σ²_∞(x) +
  (1/N) · Σ_{h: max_p_h > Q} (ε_Q(h) · li_2(x))² where
  ε_Q(h) = |S_∞(h) − S_Q(h)|. Empirically verified at 16 (anchor, x_j,
  Q) cells, predictions within 5–25% of σ_emp.
- **Knee scaling Q* ≈ √x / log x.** Empirical at 10⁷: knee_Q = 200,
  knee_max_p = 199 (predicted √x/log x = 196). At 10⁸: knee_Q ∈
  {1000, 2000}, knee_max_p ∈ {599, 1009} (predicted = 543).
- **Sharp half-Gaussian shape transition Q=100 → Q=200.** Cross-h KS
  median p-value jumps from 0.0015 (Q=100, 0/10 cells pass p>0.1) to
  0.96 (Q=200, 10/10 cells pass).
- **Cost dimension.** S_Q(h) by trial division costs min(Q, √h) odd
  integers tested per h. Average cost across 26-h ensemble: 7.0 (Q=30)
  → 8.5 (Q=∞), only 18% saving. Q-truncation provides no asymptotic
  cost saving in the original P2 regime (h ≤ poly log x ⇒ √h ≤ poly
  log x ≤ Q* = √x/log x).

**Structural conclusion (S420):** Q-truncation is descriptive, not
algorithmically novel — the named-exponent decomposition cleanly
separates intrinsic noise from truncation contamination, and the
knee Q* ≈ √x/log x is *algebraic* in x, not polylog. The original
P2 polylog regime (h ≤ poly log x) bypasses Q-truncation entirely:
S_∞(h) computation is √h ≤ poly log x. **Polylog HL evaluator with
named-exponent error ε(x, h-ensemble) ≈ F · √(S̄ · li_2(x)) (matching
slot-2's empirical 0.36–0.70 F-factor) is the slot-4 wrap target.**

**Slot 4 (S421) result — conditional theorem T-A' + Corollary B'.**
Theoretical wrap. The slot-3 empirical variance decomposition is lifted
to a CONDITIONAL THEOREM under the cross-h Hardy–Littlewood Random-
Residual Hypothesis (HLRH). HLRH consists of three asymptotic
statements on cross-h ensembles at fixed x: (a) first-moment vanishing,
(b) second-moment asymptotic to F²_H · S̄_H · li_2(x), (c) cross-h
decoherence. HLRH is the cross-h analogue of the S195 random-phase
model and the S244 Montgomery pair-correlation hypothesis transposed
from K-axis (zero index) to h-axis (gap index). It is implied by the
strong k-tuple HL conjecture for k = 4 plus GUE-side correlations on
prime-pair-count joint distributions (currently unproven beyond
Bombieri–Davenport linear-bound levels).

- **Theorem A' (T-A'):** under HLRH(H), the cross-h sample variance
  σ²_HL_Q(x, H) = (1+o(1))·F²_H·S̄_H·li_2(x) + (1/N)·Σ_{h: max_p_h>Q}
  ε_Q(h)²·li_2(x)² + o(S̄_H·li_2(x)) for any Q ∈ ℕ ∪ {∞}.
- **Knee Corollary:** Q* ≍ √x/log x is the quadrature-equality solution
  of σ²_∞ = σ²_truncation. Empirically validated at 10⁷ (predicted
  196, knee 199–200) and 10⁸ (predicted 543, knee 599–1009).
- **Corollary B' (algorithmic):** for h-batches H with max h ≤
  (log x)^β and |H| ≤ poly(log x), the full HL evaluator costs
  O(N·√(max h)) = poly(log x) per batch and has cross-h L²-typical
  error ε_typ(x, H) ≤ (1+o(1))·F_H·√(S̄_H·li_2(x)) ≍ F_H·√(S̄_H·x)/
  log x — the named-exponent error bound.
- **K-axis ↔ h-axis analogy table:** §8 of slot4_theorem.md provides
  the structural correspondence with S244's proof; the bilinear-form
  machinery is the cross-h transposition of Goldston–Montgomery 1987.

**Structural conclusion (S421):** Thread 8 closes
DONE_PARTIAL_POSITIVE_CONDITIONAL. The aggregate Thread 8 contribution
is a polylog-time HL evaluator with named-exponent cross-h L²-typical
error, conditional on HLRH. This is the project's **second** A-shape
positive-direction CONDITIONAL theorem on an adjacent π-related
computation (after Thread 7 / S244), and the **first** on the h-axis
(gap family). The Thread 8 wrap is structurally parallel to Thread 7's
wrap but weaker in two specifics: (i) F_H is ensemble-dependent rather
than universal (lacking Thread 7's flat decade-stability), (ii) HLRH
is a stronger / less-studied hypothesis than RH + Montgomery.

**Files:** `experiments/analytic/batched_pi_h/slot4_theorem.md`,
`experiments/analytic/batched_pi_h/slot3_q_truncation_results.md`, and
`archive/sessions/session421_commit_p2_theorem_wrap.md`.

### P3 — Approximate π(x) with bounded error in polylog time [Thread 7 ACTIVE]
**Problem:** what is the smallest ε(x) such that there is a polylog
algorithm computing π(x) ± ε(x)? The Riemann smooth approximation
R(x) gives ε ≈ √x. Can polylog buy ε = O(x^{1/2-δ}) for some δ > 0?

**Status: CLOSED-CONDITIONAL — Thread 7 done at S244, all 5 slots
complete. Aggregate result: polylog-time algorithm for approximate
π(x) with named-exponent error ε(x) ≤ √x · log log x / log^β x for
any β > 1, **conditional on RH + Montgomery's pair-correlation
conjecture**. Empirically verified across 3 decades (S241), kernel-
optimal across 17 kernel families (S242 + S243, 180 cells), rigorised
modulo Montgomery (S244). See `experiments/analytic/polylog_approx_pi/
slot5_theorem.md` for the conditional theorem statement and proof
outline.**

**Why partial-positive plausible:** R(x) + K zeros gives, under
Montgomery random-phase heuristic (S195 variance formula),

  σ(K, x) ≈ √x · log K / (π √(2K) · log x).

For K = (log x)^α the named-exponent corollary (S240) is

  σ(x, K = log^α x) ≈ α · √x · log log x / (π √2 · log^{1 + α/2} x).

Inverting: for any β > 1, taking K = log^{2(β-1)} x zeros gives a
polylog-time algorithm with typical error
ε(x) ≤ √x · log log x / log^β x. (Note: this CORRECTS an earlier
version of this entry that claimed K = log²x gives
ε ≈ √x · log log x / log⁴ x — that used K rather than √K in the
denominator; correct formula gives ε ≈ √x · log log x / log² x at
K = log²x.)

**S240 empirical confirmation** (slot 1, single-anchor measurements
x ∈ {10⁵, 10⁶, 10⁷, 10⁸, 10¹⁰}, K up to 8000 zeros): median
empirical / σ-predicted ratio = 0.476, mean 0.554 — consistent with
half-Gaussian factor √(2/π) ≈ 0.798 modulated by GUE pair-correlation
reduction (~0.74 from S195). At x = 10¹⁰, K = 8000 (≈ log³x):
|err| = 48 vs √x = 100000, i.e., ε/√x ≈ 5×10⁻⁴.

**S241 empirical confirmation** (slot 2, multi-sample N=30 per anchor
at x ∈ {10⁷, 10⁸, 10⁹}, 360 (x, policy) data points): half-Gaussian
**shape** of |err|/σ_eff confirmed — median KS p-value 0.69 across 9
(anchor, K ≥ log²x) cells. The GUE pair-correlation factor
σ_eff/σ_pred = 0.755 ± 0.06 stable across 3 decades, extending S195's
0.74 figure to x = 10⁹. ε/√x at K=8000 fixed: 9.7e-4 → 8.8e-4 →
7.8e-4 from 10⁷ to 10⁹, cross-decade ratio 1.25 vs predicted 1.286
(1/log x scaling, agreement within 3%). Theoretical extrapolation
table extended to x = 10²⁴: at x = 10¹⁵ with K = log⁴x ≈ 9.3×10⁵
zeros, σ/√x ≈ 7.7×10⁻⁵.

**S242 slot-3 negative-shape closure** (smoothed kernel selection):
N=20 paired samples × 9 compactly-supported kernels (Hann, Hamming,
Triangle, Riesz, Riesz⁴, Tukey-25, Tukey-50, Cosine, Hard) × 4 K-values
× 3 anchors → **0 of 96 (anchor, K, kernel) cells show smoothed kernel
beating hard cutoff at p < 0.05**. Mean σ_eff(kernel)/σ_eff(hard)
across 12 (anchor, K) cells, by kernel: tukey25=1.04, tukey50=1.11,
cosine=1.14, riesz=1.12, triangle=1.23, riesz4=1.21, hamming=1.20,
hann=1.23. Decisive hard wins at 10⁸/K=2000+4000 and 10⁹/K=2000+4000:
ratios 1.06-1.62, paired sign test 12-17/20 hard wins, p < 0.01 in
4 cells. **Closes S202 wrap §"Non-Gaussian smoothing kernels"
legitimate-falsifier**: hard cutoff is empirically L2-optimal for the
partial-sum approximation π_K(x) under the GUE-corrected random-phase
regime. The slot-1 named-exponent corollary σ ≈ √x · log K /
(π√(2K) · log x) is **kernel-optimal** in the compactly-supported
family.

**S243 slot-4 negative-shape closure** (paired / non-symmetric /
position-correlated kernels): N=20 paired samples × 7 paired families
+ hard baseline × 4 K-values × 3 anchors → **0 of 84 (anchor, K,
kernel) cells show paired-kernel beating hard cutoff at p < 0.05**.
Smallest p_kernel_beats_hard = 0.252 (boundary_pair, 10⁷ K=500). Mean
σ_eff(kernel)/σ_eff(hard) by family: boundary_pair=0.999,
half_int=0.999 (statistically indistinguishable from hard);
paired_riesz=1.12, paired_triangle=1.23, paired_hann=1.23 (similar to
slot-3 symmetric); antipair_03=2.12, antipair_05=3.24 (catastrophically
worse — antipair_05 hits 5.05× hard at K=4000, 10⁹). The antipair
scaling Var(antipair)/Var(hard) = 1 + (1−w)² K matches the random-phase
prediction exactly. **Closes the LAST remaining kernel-axis falsifier**
not addressed by S196 (bandwidth-axis) or S242 (symmetric kernels):
**hard cutoff is empirically L2-optimal across all 17 kernel families
tested (8 symmetric + 7 paired + hard, 180 cells across slots 3+4)**.
F_GUE := σ_eff²/σ_pred² = 0.55 ± 0.06 stable across all kernels
(GUE pair-correlation is kernel-invariant at second-moment level).

**S244 slot-5 theoretical wrap (conditional theorem):**
`slot5_theorem.md` adapts Goldston–Montgomery 1987 ("Pair correlation
of zeros and primes in short intervals") bilinear-form analysis to
the truncated-zero-sum test function. **Theorem A (slot 5,
conditional):** under RH + Montgomery, for H ∈ [X^ε, X log^{−2}X]
and K ∈ [log²X, X^{1−ε}], (1/H) ∫_X^{X+H} (π(y) − π_K(y))² dy =
(1+o(1)) · X · log²K / (2π² K · log²X). **Corollary B (algorithmic;
conditional):** under same hypotheses, K = ⌈(log x)^{2(β−1)}⌉ gives
polylog-time algorithm with L²-typical error ε_typ(x) ≤ (1+o(1)) ·
(β−1) √2 · √x · log log x / (π · log^β x). The valid range (★) is
explicit. Montgomery's conjecture is used only for the *close-pair*
off-diagonal bound; the far-pair bound is RH-only. Under RH alone
the same proof gives σ²_RH ≤ X · log²K · log²log K / (2π² K · log²X)
— same exponent in log X, log²log K factor weaker.

**A-grade target (achieved as conditional):** the named-exponent
corollary is now a precise CONDITIONAL theorem under RH + Montgomery,
not just heuristic. **Unconditional version remains open** —
Montgomery's pair-correlation conjecture is itself unproven, and
under just RH the named exponent is preserved up to a (log log x)²
factor in σ.

**Outstanding open directions for future threads:** (i) prove a
Wigner-repulsion sub-conjecture sufficient for the close-pair bound
without the full Montgomery conjecture; (ii) replace the L²-typical
window-average with a pointwise-in-y bound (the half-Gaussian tail
of S241 suggests pointwise error is √(log K) larger than typical
at the worst x); (iii) push to non-linear post-processing of {c_k}
(empirical Bayesian shrinkage) — outside the slot-3+4 closed linear
framework; (iv) push to adversarial per-x K-policy selection
(requires a second-moment estimator at fixed x).

**Budget:** 5-session commit thread (Thread 7). DONE — all 5 slots
complete at S240, S241, S242, S243, S244.

### P4 — Twin-prime / k-tuple count in narrow windows [Thread 9 CLOSED-CONDITIONAL_PARTIAL_POSITIVE at S433; all 5 slots done]
**Problem:** `π_H(x; w) = #{p ∈ [x, x+w] : p+h prime for all h ∈ H}`
for admissible H = {0, h_1, ..., h_{k-1}} ⊂ ℕ and narrow window
w = polylog(x). Batched on x_i.

**Status: CLOSED-CONDITIONAL_PARTIAL_POSITIVE — Thread 9 done at S433.**
Slots 1-5 complete:
- **Slot 1 (S429, B):** sieve-shared batched-x evaluator; 5-8× speedup
  at M=64.
- **Slot 2 (S430, B):** F_HL_kt = 0.87 ± 0.03 across 3 decades.
- **Slot 3 (S431, B):** pair-correlation derivation; F_pred matches
  F_emp <0.5% at 10⁸ wide cell; Theorem 1 (Gallagher-Poisson under HL).
- **Slot 4 (S432, B):** asymptotic shape correction; structural
  candidate Δ ∼ -12 C₂² · w · log w · log log w (matched empirical
  -5.36 to 2.4% on K up to 200K).
- **Slot 5 (S433, B, FINAL WRAP):** REFUTED slot 4's -12 C₂² as a
  unique fit on extended K=350K data — Model B coefficient varies
  -5.5 to -8.9 across reasonable fit ranges (slot 4's -5.36 was a
  16-cell local minimum); first exact full single-prime decomposition
  S_1 captures 25-29%; first explicit cross-prime second-order S_2
  (primes ≤ 200) captures 16-35%; combined 41-61%, remainder 39-58%
  from cross-prime tail + higher-order. Thread 9 partial-positive
  characterization with conditional Theorem 1 from slot 3 +
  empirical leading coefficient interval [5.0, 9.0] (refined from
  slot 4's point estimate -12 C₂²).

**Conditional theorem (Thread 9 wrap):** under HL-4 + slot 5
Conjecture (Δ(w) = -A_∞(w)·w·log w·log log w·(1+o(1)) with
A_∞ → A_* ∈ [5.0, 9.0]):
```
F²(x; w) = Var[N_2(x; w)] / E[N_2(x; w)]
         = 1 + Δ(w) / (C_2 · w · log²x) + ε_GM(x, w)
```
where ε_GM is the unaccounted Goldston-Montgomery zero contribution
(slot 3/4 systematic ±0.003 residual). For w/x → 0 with x → ∞:
F²(x; w) → 1.

**Six falsifiers** documented in slot 5 wrap (slot5_thread9_wrap_results.md).

**Successor proposals (PROPOSED):**
- **P4-extension-a:** Goldston-Yıldırım 2007 partial-sum derivation
  on HL 4-tuple singular series. If derived rigorously, resolves
  slot 5's open structural question.
- **P4-extension-b:** GM zero-residual analysis at slot 3 cells with
  slot 4 corrected evaluator (re-run F_pred at 10⁶ wide, 10⁷ wide;
  test ε_GM ~ const · log w / log x at 3 cells).

Slot 1 (S429): sieve-shared batched-x k-tuple primitive verified correct
at 72/72 cells; speedup ratio T_batched/T_naive at M = 64 measured at
0.188-0.213 at x = 10⁶ (5×) and 0.119-0.135 at x = 10⁷ (8×) for
correlated narrow-window distributions, k-independent across k ∈ {2
twin, 3, 4 admissible} within ±10%. Uncorrelated distribution gives
9-49× WORSE than naive (anti-amortising shared-sieve over Θ(x_max/2)
range — same structural shape as Thread 5). HL approximation
HL_H(x; w) = C(H)·w/log^k x matches empirical at 6 cells within
≤ 0.34σ_Pois.

Slot 2 (S430): cross-x HL residual distribution at fixed admissible H,
N=200 disjoint windows per cell, 18 cells (3 anchors × 3 H × 2 w-regimes).
**F_HL_kt = σ_eff/σ_pois = 0.87 ± 0.03 for k=2 across 3 decades**
(range 0.09 wide / 0.045 narrow), MORE decade-stable than Thread 8's
F_HL ∈ [0.36, 0.70] (range 0.34) by ~5×. Half-Gaussian KS passes at
k=2 wide regime: 2/3 anchors p > 0.1, 3/3 p > 0.05. For k≥3 wide:
F → 1.0 (Poisson-exact at variance level, range 0.025 across 3 decades).
**Methodological correction to S419**: cross-x at fixed H IS iid-like
for *windowed* counts, where Thread 8 cross-x at fixed h FAILED for
*cumulative* π_h.

**Why partial-positive plausible:** narrow-window k-tuple count under
the Hardy-Littlewood conjecture is `~ C(H) · w / log^k x`. Per-x
amortisation possible across correlated narrow-window x_i via shared
segmented sieve.

**Slot 1 (S429) result — dichotomy transposes, magnitude smaller.**
Thread-5 / S224 correlation dichotomy direction is preserved
(correlated → speedup grows with M and x; uncorrelated → no
shared-sieve help). Magnitude is smaller (5-8× vs Thread 5's 33×) because
per-query baseline is segmented-sieve √x · log log x · w (not Lucy DP
O(x^{2/3})). Predicted-and-confirmed √x · log log x scaling across the
10⁶ → 10⁷ measurement.

**Concrete first step:** DONE at S429. See files
`experiments/analytic/k_tuple_batched/{slot1_baseline.py,
slot1_baseline.csv (72 rows), slot1_hl_compare.csv (6 rows),
slot1_run.log, slot1_baseline_results.md}` and
`archive/sessions/session429_commit_p4_baseline.md`.

**A-grade target:** named-exponent error bound for HL_H(x; w) on the
cross-x ensemble at fixed admissible H, conditional on a precise
HLRH-x hypothesis — analogue of Thread 7 Corollary B (K-axis) and
Thread 8 Corollary B' (h-axis), transposed to the cross-x ⊗ k-tuple
axis. The empirical measurable speedup over per-window sieving is
already established at slot 1 (5-8× growing with √x).

**Slot 3 (S431) result — pair-correlation derivation matches slot 2.**
Verified the prime-by-prime cancellation identity ⟨S_4 factor at p⟩_uniform_m
= S_2² factor at p exactly (max deviation 4×10⁻¹⁶ for first 168 primes).
Combined with admissibility factors 2 (p=2) and 3 (p=3), this yields
⟨S_4⟩_admissible = 24 C_2² ≈ 10.4596. Theorem 1 (windowed-twin
Gallagher-Poisson under HL): Var[N_2(x;w)]/E[N_2] → 1 as x→∞ with
w/x→0. Direct computation of T(w) for slot-2 cells gives sub-leading
Δ(w) = T(w) - 2C_2²w² ≈ -5.72·w·log(w) + 24C_2²·w (intercept matches
24C_2² to 0.35%). F²_pred(x;w) = 1 + Δ/(C_2·w·log²x) matches slot 2
F²_emp to **0.2% at x=10⁸ wide**, **1.3% at 10⁷ wide**, **3.4% at 10⁶
wide**; monotone improvement with x. Narrow-regime cells deviate
4-7% (discrete-count quantization). Residual F_pred − F_emp > 0 in all
6 cells indicates an unaccounted ~0.002-0.07 Goldston-Montgomery
zeros contribution beyond HL singular series.

**Slot 4 (S432) result — asymptotic shape correction + slot 3 bug.**
Extended T(w) computation to K = 200,000 (w = 1.2 × 10⁶), 22 cells
across log w ∈ [5.7, 14.0]. Found:

(i) **Slot 3 software bug** (off-by-one in tail handling): skipped
    the first prime > diam(H) in singular-series evaluation, biasing
    S_4 high by 1/(1−4/p_tail). Fast slot 4 evaluator validated to
    bit precision against Hardy-Littlewood prime quadruplet constant
    S_4(0,2,6,8) = 4.1511 (slot 3 had 4.4571). Slot 3 cells |Δ|
    underestimated by 4–6%.

(ii) **Slot 3 α ≈ 5.72 REFUTED as finite-w artefact.** Rolling-band
     α(log w) grows monotonically: 6.43 (K=50..200) → 9.73
     (K=50000..200000). Linear-in-log-w model rejected.

(iii) **Best 3-param fit** Δ/w ≈ −5.36 log(w) log log(w) + 9.30 log(w)
      − 22.37 (RMS rel err 1.5%); log² model rejected (coefficient
      flips sign as fit range moves to larger K).

(iv) **Structural candidate Δ(w) ∼ −12 C_2² · w · log(w) · log log(w)**
     with −12 C_2² = −5.230 matching empirical −5.36 to 2.4%. Single-
     prime Mertens-style heuristic captures 32% of magnitude;
     remaining 68% from cross-prime + boundary primes (analytic
     derivation incomplete).

(v) **Slot 3 cells corrected**: F_pred at 10⁸ wide cell shifts from
    0.9137 (slot 3) to 0.9080 (slot 4 corrected) vs F_emp = 0.9113;
    |F_diff| stays <0.5% but flips sign — overshoot becomes
    undershoot. Unaccounted ±0.003 plausibly Goldston-Montgomery
    zeros contribution.

See `experiments/analytic/k_tuple_batched/{slot4_alpha_derivation.py,
slot4_h_residual.csv, slot4_t_delta.csv, slot4_alpha_fits.csv,
slot4_slot3_comparison.csv, slot4_run.log, slot4_alpha_derivation_results.md}`
and `archive/sessions/session432_commit_p4_alpha_derivation.md`.

**Slot 5 (FINAL, S433) result:** option (a) (Goldston-Yıldırım
derivation) UNREACHABLE in slot 5 scope (would require Selberg-
Delange / GY 2007 machinery beyond current toolkit); on extended
K=350K data, slot 4's structural identification with -12 C₂² is
REFUTED as a unique fit (Model B coefficient varies -5.5 to -8.9).
Option (b) (GM zero-residual) UNDERTESTED with only 1 corrected
F_pred cell. Option (c) (Thread 9 wrap) is realised. See
`experiments/analytic/k_tuple_batched/{slot5_thread9_wrap.py,
slot5_extended_t.csv, slot5_decomposition.csv, slot5_gm_residual.csv,
slot5_run.log, slot5_thread9_wrap_results.md}` and
`archive/sessions/session433_commit_p4_thread9_wrap.md`.

**Budget:** 5-session arc DONE. Slots 1-5 all complete per
`.commit_state thread:NONE_AWAITING_USER_ESCALATION sessions_used:0
prev_thread_9:p4_k_tuple_narrow_window_batched_on_x_DONE_PARTIAL_POSITIVE_CONDITIONAL`.

---

## §2. Conditional / Heuristic Variants

### P5 — π(x) under GRH with explicit constants [CLOSED-B at S436 — Thread 10 wrapped B-NEGATIVE]
**Problem:** Galway 2004 gives `K = O(x^{1/2+ε})` zeros under GRH.
What's the explicit constant? At what x does GRH-conditional become
faster than HKM unconditional in *practice*?

**Status: CLOSED-B at S436. Thread 10 wraps as B-NEGATIVE structural
insight — the P5 tightening goal is unachievable not because of compute
limits but because the Galway-shape `K = c · √x · log²x` with c=const is
asymptotically loose. The worst-case-of-N=30 asymptotic is Thread-7-shape
`K ~ x · log²K / log²x`, validated via three independent angles across
log10 x ∈ [4.0, 5.5]: slot 1+2 direct measurement (c_emp drift 0.13 →
0.49, refuting const c), slot 2 extended K_max=20000 measurement at log10
x = 5.5 (rigorous LB c_emp > 0.222 + worst |err|@K=20000 = 1.609 + err ~
1/√K extrapolation giving c_emp ≈ 0.574), and slot 3 RMS-based σ_pred
cross-check (σ_eff(K=20000)/σ_pred(K=20000) = 0.74-1.05 across 3 high-x
anchors, mean 0.897 ± 0.16 — independent of slot 2's worst-case extrapolation).
GUE-corrected Thread-7-shape predicts K_emp(5.5) = 44,341, K_emp(6.0) =
145,436. Literature audit across {Lagarias-Odlyzko 1987, Galway 2004,
FKBJ 2017, Büthe 2018, Platt-Trudgian 2012} confirms no published source
gives explicit small numerical c for the unsmoothed Riemann-R partial
sum K-budget; slot 1+2+3 is the FIRST empirical prefactor measurement
and FIRST cross-decade refutation of const-c at the worst-case tail.

**Slot 3 (S436) added Parts A, B, C:**

| Part | Content | Headline result |
|---|---|---|
| A | RMS-based σ_eff/σ_pred at K=20000 | mean 0.897 ± 0.16 across 3 anchors — Thread-7 σ_pred validated |
| B | GUE-corrected Thread-7 K_emp predictions | K_emp(5.5)=44341, K_emp(6.0)=145436 |
| C | Empirical c_emp vs Thread-7 (13 anchors) | log10 x ≤ 4.6 within 17%; log10 x ≥ 5.0 K_max-limited |

See `experiments/analytic/galway_constant/{slot3_thread7_validation_results.md,
slot3_literature_audit.md}` and `archive/sessions/session436_commit_p5_galway_thread10_wrap.md`.

Slot 1 measured empirical `c_emp` in `K = c · √x · log²x` at log10 x ∈
{4.0, 4.5, 5.0} (K_max=8000): 0.177 / 0.210 / >0.191 LB. Slot 2
extended the measurement on two paths:

**Slot 2 path (a) — finer x-grid (K_max=8000, 10 anchors):** c_emp
mean 0.151 ± 0.044 across log10 x ∈ [4.0, 4.9], range [0.106, 0.249].
Single-decade dynamic range insufficient to distinguish Galway-shape
from Thread-7-shape within sample-variance noise. See
`experiments/analytic/galway_constant/slot2_finegrid_results.md`.

**Slot 2 path (b) — extended K_max=20000 at log10 x ∈ {5.0, 5.3, 5.5}:**
combined zeros database with 12,000 new zeros (k=8001..20000 via 12
parallel mpmath.zetazero workers).

| log10 x | x | K_emp | c_emp | c_emp_T7_pred |
|---:|---:|---:|---:|---:|
| 4.00 | 10000 | 1250 | 0.1474 | 0.1785 (slot 2 finegrid) |
| 4.50 | 31623 | 4750 | 0.2488 | 0.2638 (slot 2 finegrid) |
| 5.00 | 100000 | 9000 | 0.2147 | 0.3942 (slot 2 ext) |
| 5.30 | 199526 | 9000 | 0.1353 | 0.5044 (slot 2 ext) |
| 5.50 | 316228 | **>20001** | **(LB > 0.222)** | **0.5958** (slot 2 ext) |

**Headline slot 2 result: Galway-shape c_emp = const ≤ 0.222 REFUTED
at log10 x = 5.5.** K=20000 budget is INSUFFICIENT for ε=1 worst-of-30
even at the upper end of the empirical c_emp range. Worst sample's
|err| at K=20000 is 1.609; extrapolation under err ~ 1/√K gives K_emp
≈ 51,778 → c_emp ≈ 0.574 — matches Thread-7-shape prediction
c_emp_T7(5.5) = 0.596 within 4%. Slot 2 σ_eff/σ_pred ratio at K=20000
across 3 anchors: 1.0 / 0.83 / 1.0 (mean 0.95), extending Thread 7's
ratio-drift trajectory: 0.72 (10⁴) → 0.79 (10⁴·⁵) → 0.88 (10⁵) → 1.0
(10⁵·⁵). See `experiments/analytic/galway_constant/slot2_extended_results.md`.

**Implication: E6.1 (Galway 2004 K = O(x^{1/2+ε})) is asymptotically
loose at the worst-case-of-N tail.** Actual K_emp scales like Thread-7-
shape `K ~ x · log²K / log²x`, not `K ~ √x · log²x`. The "c ≈ 0.20"
constant from slot 1 is a finite-x phenomenon at log10 x ∈ [4, 5], not
asymptotic. **The slot 2 path-(b) measurement at log10 x = 5.5 is the
project's first cross-decade refutation of a Galway-style worst-case
explicit-formula prefactor at the worst-case-of-N tail.**

**Why partial-positive plausible:** the constant has never been
empirically tightened. Slot 1+2 establish that the constant is
not actually asymptotically constant — Thread-7-shape dominates.

**A-grade target (NOT achieved as positive):** GRH-conditional π(x)
algorithm strictly faster than HKM on a real benchmark. Slot 1+2
empirical `c_emp ≈ 0.16-0.25 at log10 x ∈ [4, 5]`; whether this beats
HKM/primecount in practice depends on per-zero cost (R_at_rho ~0.43ms
at dps=18 vs HKM ~1µs/element) — NO at small x. The Thread-7-shape
asymptotic super-Galway scaling further weakens any practical
advantage at large x.

**Slot 3 (S436) DONE — Thread 10 wrapped B-NEGATIVE:**
- Path (a) executed: literature audit + RMS-based Thread-7 cross-check.
- Slot 3 Part A: σ_eff(K=20000)/σ_pred(K=20000) ratio = 1.05 / 0.74 / 0.90
  across log10 x ∈ {5.0, 5.3, 5.5}; mean 0.897 ± 0.16 — independent
  empirical validation of Thread 7's σ_pred at log10 x ≥ 5.0.
- Slot 3 Part B: GUE-corrected Thread-7 K_emp(5.5) = 44341 (2.2× the
  K=20000 budget; consistent with slot 2's |err|@K=20000 = 1.609).
- Slot 3 Part C: empirical/Thread-7(f=0.755) ratio = 1.01-1.32 at log10
  x ∈ [4.0, 4.5]; K_max-limited at log10 x ≥ 5.0.
- Literature: no published source gives explicit small numerical c for
  the unsmoothed Riemann-R partial sum K-budget; slot 1+2+3 is FIRST
  empirical prefactor measurement and FIRST cross-decade const-c
  refutation. Path (b) NOT executed (would convert extrapolation to
  direct measurement; not needed for closure).

**Final outcome:** B-NEGATIVE structural insight (P5 tightening goal
unachievable) + B-PARTIAL-POSITIVE (Thread-7-shape worst-case-of-N
asymptotic established across 1.5 decades).

**Budget:** 3-session arc. Slots 1+2+3 done at S434+S435+S436. Thread 10
CLOSED.

### P6 — Conditional polylog under stronger heuristics
**Problem:** under Hardy-Littlewood + Cramér + Montgomery joint
heuristics, is polylog π(x) achievable for *typical* x (i.e., on a
density-1 set, with a measurable exception set)?

**Why partial-positive plausible:** the probabilistic-prime-model
uncertainty around π(x) might be polylog-extractable for "average"
x even if worst-case is not.

**Concrete first step:** Cramér-style sampling argument. Pick x at
random from [10⁹, 10¹⁰], use only K = log³ x zeros, measure failure
rate.

**A-grade target:** conditional polylog π(x) on a density-1 subset of
x with explicit exception bound.

**Budget:** 3-session arc.

---

## §3. Restricted / Sparse Variants

### P7 — π(x) for x of restricted form [CLOSED-B at S246, per-query side]
**Problem:** for x = 2^n or x = 10^n or other "structured" x, can π(x)
be computed in polylog?

**Status: per-query side CLOSED-B at S246.** Three independent
structural tests on `(π(2^k))_{k=1..56}` (OEIS A007053) against Monte
Carlo random nulls all return null-baseline values:
- BM linear complexity of `sign(π(2^k) − R(2^k))` over GF(2) = 28
  vs MC mean 28.25 ± 1.01.
- Max lag-1..10 autocorrelation of `(π(2^k) − R(2^k)) / 2^{k/2}`
  = 0.283 vs MC p999 = 0.437 (single-lag p99 = 0.328; not significant
  after 10-lag Bonferroni).
- sympy.primepi (Meissel-Lehmer) cold-start subprocess timing ratio
  `T(2^k) / T(2^k ± 1)` = 1.013 averaged over k ∈ {28, 30, 32}.

**Why P7 closes negatively** (per-query side): at x = 2^k the
explicit-formula phases `γ_n · k log 2 mod 2π` are Weyl-equidistributed
in k for every zero `γ_n` independently of the dyadic form, and the
Lucy / Meissel-Lehmer outer loop processes `{n ≤ √x}` and
`x^{2/3}-smooth` counts which carry no binary-representation side-
channel. Per-query dyadic structure is genuinely orthogonal to the
HKM cost.

**Budget spent:** 1 session (S246), of 2 estimated. Per-query side
fully closed; batched-on-k side moves to a successor entry.

**See:** `experiments/wildcard/dyadic_pi_structure/`,
`archive/sessions/session246_f6_dyadic_pi_structure.md`,
NOVELTY_CHALLENGES.md §F6 (which now hosts successor challenges).

### P7.b — Batched dyadic π queries (cross-k amortisation, NEW S246)
**Problem:** for `x_1 = 2^{k_1}, ..., x_M = 2^{k_M}` with shared zero
database `(γ_n)_{n=1..K}`, can the per-query amortised cost over the
M dyadic anchors be polylog?

**Why partial-positive plausible:** the zero table is k-independent.
Computing all `Re(R(x_i^ρ_n))` for i = 1..M and n = 1..K costs
O(M·K) total (or O(K log M) if exploiting cosine sums via NUFFT).
Per-query amortised O(K). For K = polylog(max x_i), per-query =
polylog. **This is exactly Thread 5 / S224 Correlation Dichotomy
transposed to the dyadic family.**

**Concrete first step:** evaluate `R_K(2^{k_i})` for k_i = 1..56 with
shared zero table at K ∈ {log²x, log³x}, measure (i) per-query
amortised wall-clock and (ii) accuracy versus exact π(2^{k_i}).
Compare to Thread 7 (S240-S244) error scaling.

**A-grade target:** measurable per-query speedup over per-x HKM at
M ≥ poly(log max x).

**Budget:** 2-session arc.

### P8 — Sparse-precision π queries (batched on precision, not x)
**Problem:** compute the first k bits of π(x) for varying k. Single-k
is fixed cost; batched on increasing k might amortise.

**Why partial-positive plausible:** R(x) at low-precision gives early
bits cheaply; refining requires more zeros. Per-bit amortised cost
might be sub-linear in total bits.

**Concrete first step:** profile per-bit Galway evaluator at x = 10⁹
for k ∈ {10, 20, 40, 80, 160} bits.

**A-grade target:** strictly sub-linear per-bit cost.

**Budget:** 2-session arc.

---

## §4. Side-channel / Quantum Variants

### P9 — Quantum batched π
**Problem:** Thread 5 closed *classical* batched amortisation.
Quantum batched π queries might admit exponentially-faster amortised
cost if zeros can be encoded in quantum state.

**Why partial-positive plausible:** Shor's algorithm pattern —
quantum can sometimes give exponential amortisation that classical
can't. CTQW closure (E7.20) was for single-x; batched is open.

**Concrete first step:** literature survey of quantum-NT primitives
for π(x). Identify whether any quantum amortisation argument applies.

**A-grade target:** quantum batch-polylog π(x) algorithm.

**Budget:** 3-session arc (heavy literature, possibly inconclusive).

### P10 — Adaptive π queries
**Problem:** π(x_1) computed first, then x_2 chosen adaptively as a
function of π(x_1). Aggarwal binary-search is the simplest case.
What about non-binary-search adaptive strategies?

**Why partial-positive plausible:** Thread 5 closed Aggarwal as
suboptimal-amortisable. Other adaptive strategies (adversarial,
information-theoretic-optimal) might amortise better.

**Concrete first step:** information-theoretic lower bound on
adaptive π-query complexity. Compare to Aggarwal.

**A-grade target:** adaptive query strategy strictly better than
Aggarwal on a concrete benchmark.

**Budget:** 3-session arc.

---

## §5. Incidence-Geometric Variants

### P11 — Minimum number of straight lines covering all primes ≤ N under a 2D embedding [Thread 11 ACTIVE]

**Problem:** under a chosen 2D embedding `Φ: ℕ → ℝ²` (Ulam spiral,
residue-class grid, polynomial-image grid), let `L_Φ(N)` = the minimum
number of straight lines required to cover all points `Φ(p)` for
primes `p ≤ N`. Trivially `L_id(N) = 1` for the identity embedding on
ℤ (primes are collinear). What is `L_Φ(N)` for non-trivial Φ?

**Why partial-positive plausible:**
- The **Ulam spiral** (Stanisław Ulam, 1963) plots integers along a
  square spiral and visually shows *diagonal lines of primes*. Each
  such diagonal line corresponds to a quadratic form `f(n) = an² + bn
  + c` that produces unusually many primes. Hardy-Littlewood's
  prime-tuple conjecture for quadratics gives explicit asymptotic
  density.
- If `L_Ulam(N)` admits a named-exponent scaling `L_Ulam(N) ~ π(N)^α`
  for some `α < 1`, that's a non-trivial *incidence-geometric prime-
  density theorem* — the project's first such result.
- Known incidence-theoretic baseline: by **Szemerédi-Trotter (1983)**,
  for `n` random points and `m` lines, the number of incidences is
  `O((nm)^{2/3} + n + m)`. For random points in [1, N]², covering
  needs Ω(n^{1/2}) lines; primes might do strictly better via HL
  alignment.

**Cross-domain ingredients:** Szemerédi-Trotter incidence theorem
(1983 *Combinatorica* 3); Erdős extremal incidence problems; Ulam
spiral (Ulam 1963 *Math. Mag.* 36); Hardy-Littlewood prime-tuple
conjecture for quadratic forms (Conjecture F, 1923); Stanley 1989
*Adv. Math.* (matroid-theoretic line cover lower bounds).

**A-grade outcomes:**
- (a) **Named-exponent `L_Φ(N) ~ π(N)^α` with α < 1** for one or more
  embeddings, with HL-backing for the dominant lines → first explicit
  incidence-geometric theorem on prime structure
- (b) **Polylog-time approximation of `L_Φ(N)`** via line sampling
  exploiting HL structure on quadratic forms — concrete algorithmic
  speedup
- (c) **Matched-baseline z-score ≥ 5σ** showing primes are
  non-trivially compressible in lines compared to random points of
  same density on the embedding

**Pre-stated falsifiers:**
- (E mode B-grade) `L_Ulam(N) / L_random(N) → 1` as N grows; HL
  alignment provides only constant-factor compression. This is the
  predicted closure consistent with the project's wider HL-everywhere
  pattern.
- (I mode B-grade) `L_Ulam(N)` scales like `π(N)^{2/3}` (random
  Szemerédi-Trotter lower bound), no HL gain.
- (A mode) any deviation from (E) or (I) above with structural
  reason — would constitute the partial-positive A-grade outcome.

**Concrete first step (slot 1):** brute-force compute `L_Ulam(N)` at
N = 10⁴, 10⁵, 10⁶ via greedy line-detection or LP relaxation (Hough
transform variant). Compare to matched-baseline (random points on
spiral with density `π(N)/N`). Compare top-K dominant lines to
Hardy-Littlewood quadratic predictions (Euler's `n²+n+41` line is
the canonical example).

**Slot plan:** see `.commit_state` thread_11_slot_plan.

**Predicted failure mode:** likely E (lines reduce to HL-conjectured
quadratic-prime sequences, which are computable but algorithmically
already characterised — same wall, different shape). 10-15% A-grade
probability if HL alignment provides asymptotic compression beyond
random-point baseline.

**Budget:** 5-session commit thread.

---

## How to use this file

- When `commit` mode is invoked without an active thread, agents pick
  from §1 first (most plausible partial-positives), then §2-§4.
- When `reframe` mode fires (after a recent negative-shape closure),
  agents look for an entry in this file that's an adjacent variant.
- When a P_x target produces an A-grade result, mark it CLOSED-A and
  add a reference. When it closes at B (no positive), mark CLOSED-B
  and document why the predicted positive didn't materialise.
- When a session identifies a NEW positive target not on this list,
  append it under the appropriate §.

The framework's job, going forward, is to maximise the rate of
A-grade outputs from this list. The Correlation Dichotomy is the
prototype; produce more like it.

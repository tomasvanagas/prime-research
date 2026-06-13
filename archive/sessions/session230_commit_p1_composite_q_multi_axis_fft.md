# Session 230 — commit thread 6 slot 4: composite-q multi-axis FFT primitive + end-to-end zero finder

**Date:** 2026-04-29
**Mode:** commit (Thread 6 / π in arithmetic progressions, batched on
modulus q, slot 4 of 5)
**Self-grade:** **B** — substantive empirical refinement and architectural
extension of slot 3 to composite squarefree odd q via CRT-based
multi-axis FFT primitive. Multi-axis correctness verified to 1e-15;
end-to-end zero finder validated against mpmath at q ∈ {15, 35} to
~0.05 abs (consistent with slot-3's 0.033 at q=11). Slot-3 structural
lift preserved at composite q (FFT 1.75× DIRECT at q=1001 vs slot-3's
2.04× at q=1009); multi-axis FFT does NOT introduce asymptotic gain.
χ(-1) parity bug caught and fixed via mpmath cross-validation.

## Mission

From `.commit_state` Thread 6 slot 4 + S229 `recommended_next_action`:
**PRIORITY (b)**: composite-q multi-axis FFT (q = p₁p₂ → dual group
`Z/(p₁-1) × Z/(p₂-1)`); practical AP regime (q ∈ {15, 35, 105, 1001}).

Slot 3 (S229) built the prime-q end-to-end pipeline: FFT-shared full
AFE + Hardy Z + sign-change bracketing; cross-validated against mpmath
at q=11 to 0.03 abs; FFT 2.04× DIRECT at q=1009. **Slot-3's apparatus
relied on the CYCLIC structure of (Z/qZ)* for prime q.** Composite q
has non-cyclic (Z/qZ)* — the natural generalisation requires a
multi-axis FFT over the direct-product structure (Z/p₁Z)* × … × (Z/pₖZ)*.

Slot-4's task: build the multi-axis primitive, verify correctness,
test whether the structural lift propagates to composite q, and
measure end-to-end π(x; q, a) at the practical regime.

## What was built

`experiments/analytic/batched_q_amortisation/slot4_composite_q_multi_axis_fft.py`
(~600 lines):

### CRT decomposition layer

For squarefree odd q = p₁ · p₂ · … · p_K, the multiplicative group
decomposes via CRT:

```
(Z/qZ)*  ≅  ⊕_i (Z/p_iZ)*  ≅  ⊕_i Z/(p_i-1)
```

via primitive roots g_i mod p_i. Each n ∈ (Z/qZ)* corresponds to a
tuple (k₁, …, k_K) with n ≡ g_i^{k_i} (mod p_i). Characters factor:

```
χ_{j_1,…,j_K}(n) = ∏_i ω_i^{j_i k_i}    ω_i = exp(2πi/(p_i-1))
```

Implementation: `crt_decomp(q)` and `build_crt_log_table(q, ...)`.

### Multi-axis FFT primitive

Aggregate W indexed by tuple (k₁,…,k_K). Multi-axis DFT via
`numpy.fft.ifftn`:

```python
L[j_1,…,j_K] = phi · ifftn(W, axes=(0,…,K-1))[j_1,…,j_K]
```

phi = ∏_i (p_i - 1) = φ(q). The W tensor has shape (p₁-1, …, p_K-1,
n_t); ifftn over the first K axes produces L for ALL φ(q) chars at
once.

**Cross-method correctness** (vs direct per-character matmul):

```
q=15   (factors [3, 5],     φ=8):     err = 9.93e-16   rel = 2.47e-16
q=35   (factors [5, 7],     φ=24):    err = 2.98e-15   rel = 4.40e-16
q=105  (factors [3, 5, 7],  φ=48):    err = 2.66e-15   rel = 6.33e-16
q=1001 (factors [7,11,13],  φ=720):   err = 9.93e-15   rel = 1.13e-15
```

**Multi-axis FFT identity verified to floating-point precision.**

### Full AFE for composite q

Same two-step pattern as slot 3 but generalised:

1. Main term M(t, χ) via multi-axis FFT primitive.
2. Reflected term: for primitive χ mod q (squarefree odd q ⇒ all j_i ≠ 0),
   `Σ_{n ≤ N} χ̄(n)/n^{½-it} = conj(M(t, χ))` because χ̄(n) = conj(χ(n))
   for unitary characters of finite abelian groups and n ∈ ℝ_{>0}.
   Reflected = `W_χ · (q/π)^{-it} · ratio_gamma · conj(M)`. NO
   additional FFT — pointwise multiply.

**This identity holds for all primitive characters of finite abelian
groups, not just prime-q.** The simplification makes full AFE
essentially free given main term, just as in slot 3.

### Gauss sums via multi-axis FFT

`gauss_sum_composite(q, ...)` uses the same multi-axis FFT identity:
τ(χ_j) = phi · ifftn(e_a)[j_1,…,j_K] where e_a is indexed by
CRT-tuples. Cost ≈ phi · log(phi). Verified `|W_χ| = 1` for all
primitive χ at q ∈ {15, 35, 105, 1001}.

### Bug found and fixed: χ(-1) parity for composite q

For each odd prime p_i, ω_i^{(p_i-1)/2} = -1, so χ_i(-1) = (-1)^{j_i}.
Therefore

```
χ(-1) = ∏_i χ_i(-1) = (-1)^{Σ_i j_i}
```

**Initial implementation used the wrong formula** `(Σ_i j_i · (p_i-1)/2) mod 2`,
which only matches prime q (collapses to j₁ mod 2 since K=1). The
wrong formula gave ~1.4 abs zero error at q=15 first zero; corrected
to `(Σ_i j_i) mod 2` brings error to 0.05 — within slot-3 accuracy
range. **Caught via mpmath cross-validation**, not via slot-2/slot-3
test infrastructure (which doesn't cover composite q).

### Hardy Z + sign-change bracketing

Identical to slot 3:

```
Z_χ(t) = exp(i θ_χ(t)) · L(½+it, χ)
θ_χ_j(t) = (t/2) log(q/π) + arg Γ((1/2 + a_j + it)/2) - arg(W_χ_j)/2
```

`zero_brackets_vectorised` is the slot-3 vectorised primitive.

## mpmath cross-validation

```
q=15, χ_index=(1,1), t_max=40:
  slot-4: 2.697  5.229  8.383  10.232  11.896
  mpmath: 2.746  5.238  8.418  10.194  11.920
  abs diffs: 0.0490  0.0084  0.0351  0.0377  0.0245   max 0.049

q=35, χ_index=(1,1), t_max=30:
  slot-4: 3.593  6.774  8.206  9.218   10.860
  mpmath: 3.531  6.765  8.209  9.230   10.896
  abs diffs: 0.0621  0.0085  0.0028  0.0120  0.0356   max 0.062
```

**Slot-4 zero accuracy at composite q (~0.05 abs) is consistent with
slot-3's accuracy at prime q (0.033 abs at q=11).** Both bounded by
the slot-3 oversize-12 hard AFE truncation.

## End-to-end π(x; q, a) timing (3-run median)

| q     | factors      | φ    | K   | t_zf (ms) | t_psum (ms) | t_single (ms) |
|-------|--------------|------|-----|-----------|-------------|---------------|
| 15    | [3, 5]       | 8    | 50  | 3.80      | 0.16        | 3.96          |
| 35    | [5, 7]       | 24   | 50  | 4.90      | 0.35        | 5.25          |
| 105   | [3, 5, 7]    | 48   | 50  | 6.59      | 0.36        | 6.94          |
| 105   | [3, 5, 7]    | 48   | 200 | 37.82     | 0.56        | 38.38         |
| 1001  | [7, 11, 13]  | 720  | 200 | 152.27    | 14.31       | 166.58        |

### Headline comparison vs slot-3 prime-q baseline

```
                              slot-3 q=1009 (prime)    slot-4 q=1001 (composite)
T_AFE                         139 ms                   ~140 ms
T_partial_sum                  34.8 ms                  14.3 ms
T_total single-query          186.3 ms                 166.6 ms (89% of slot-3)
```

**11% wall-clock drop at q=1001 vs q=1009** attributable mostly to
smaller φ(1001) = 720 vs φ(1009) = 1008 (29% reduction). Per-character
cost essentially unchanged because N_AFE depends on q (not φ) and is
comparable.

## FFT vs DIRECT comparison

```
q=105,  K=50,  φ=48,  N_AFE=283:  FFT 5.09 ms vs DIRECT 4.24 ms = 0.83×
q=1001, K=200, φ=720, N_AFE=1872: FFT 145.24 ms vs DIRECT 254.39 ms = 1.75×

(slot-3 prime-q reference:
 q=101,  φ=100, FFT 0.95×;  q=1009, φ=1008, FFT 2.04×)
```

**Key empirical finding**: multi-axis FFT crossover from "loss" to
"win" around q ~ 200-500, same as slot-3 single-axis FFT for prime q.
At q=1001 vs q=1009, multi-axis is *slightly worse* than single-axis
(1.75× vs 2.04×) because numpy.fft.ifftn calls 1D FFT routines per
axis; small per-axis transforms (e.g., shape (6, 10, 12) at q=1001)
have proportionally higher constant overhead than one large 1D FFT
(size 1008 at q=1009).

**Multi-axis FFT preserves slot-3's structural lift in the composite-q
regime; does NOT introduce asymptotic improvement.**

## Edges composed / cited

- **E1.5** (per-query bit-content barrier): unchanged — slot-4's
  multi-axis primitive affects setup + AFE eval but not partial-sum
  eval, exactly per E1.5 prediction.
- **E2.1** (spectral / sieve interface): the multi-axis DFT identity
  generalises slot-2's cyclic DFT to direct-product groups via CRT.
- **E3.1** (CCM amortisation, downgraded): same setup-vs-eval
  decoupling pattern.
- **E6.6** (Aggarwal binary search): orthogonal — within-q vs cross-call.
- **E6.7** (Deléglise-Rivat per-query Θ(x^{2/3}/log²)): orthogonal pillar.
- **S224** (Correlation Dichotomy, Thread 5): slot-4 wall-clock shape
  ~1.75× constant lift at q=1001 NOT Correlation-Dichotomy-shaped (33×
  at M=64). Same shape as slot-3: "real-zeros for synthetic-cost"
  bounded constant.
- **S227** (slot 1, Thread 6): linear-in-φ(q)·K baseline.
- **S228** (slot 2, Thread 6): cyclic-DFT primitive that slot-4
  generalises to direct-product groups via CRT.
- **S229** (slot 3, Thread 6): prime-q end-to-end pipeline that slot-4
  extends to composite q with structural lift preserved.

## Cross-domain ingredient

**CRT-based multi-axis Dirichlet character DFT.** Generalisation of
slot-2's cyclic-group DFT primitive to direct-product groups
(Z/p₁Z)* ⊕ … ⊕ (Z/p_KZ)* via Chinese Remainder Theorem decomposition.
Implemented as `numpy.fft.ifftn` on multi-dimensional tensor indexed
by CRT log-tuples.

**NEW USED I** in this project. Added to CROSS_DOMAIN_TECHNIQUES.md
§10 family alongside slot-2 cyclic-DFT row.

References:
- Bober & Hiary, "Computing Dirichlet character L-values" (2017) for
  slot-2 antecedent (cyclic FFT primitive).
- Davenport, "Multiplicative Number Theory" §1 (CRT decomposition of
  (Z/qZ)*).

## Self-extension proposals

Per CLAUDE.md autonomy invariant:

**Successor 1 (slot 5 PRIORITY-a, NEW)**: Cross-conductor batches via
shared coprime tower for nested Q-family (e.g., q ∈ {15, 105, 1001}).
Investigate whether the cp_all summands can be shared structurally.
Likely closure under CRT independence but worth a careful empirical +
theoretical check; would seal Thread 6 across all four amortisation
axes.

**Successor 2 (slot 5 PRIORITY-c)**: Riemann-Siegel correction terms
for Dirichlet L (Berry 1995 / Coffey 2003). Drops N_AFE by 5-10× with
constant overhead; would re-test multi-axis FFT primitive at much
smaller AFE cost. Slot-3's S229 self-extension proposal, deferred.

**Successor 3 (slot 5 PRIORITY-d)**: Non-primitive character handling
at composite q via inducing-primitive-character zeros + Euler factor
corrections. Operationally important but orthogonal to the slot-4
cost question.

## Falsifiers (slot-4 bounded results)

- **F4.1** Multi-axis FFT does not give asymptotic improvement over
  slot-3 single-axis at composite q. Both achieve constant-factor
  lift in same crossover regime.
- **F4.2** End-to-end wall-clock dominated by AFE eval (N_AFE · n_t);
  composite q inherits same AFE cost as nearby prime q.
- **F4.3** Multi-axis FFT does not shrink AFE truncation requirement
  (balanced N = √(qt/(2π)) is q-only).
- **F4.4** Cross-conductor amortisation across Q-family NOT attempted
  in slot 4 (slot-3 leftover priority-a).

## Files written by this session

- `experiments/analytic/batched_q_amortisation/slot4_composite_q_multi_axis_fft.py` — new (~600 lines)
- `experiments/analytic/batched_q_amortisation/slot4_composite_q_multi_axis_fft_results.md` — new
- `experiments/analytic/batched_q_amortisation/slot4_composite_q_end_to_end.csv` — 5 rows
- `experiments/analytic/batched_q_amortisation/slot4_composite_fft_vs_direct.csv` — 2 rows
- `.commit_state` — sessions_used 3 → 4, session_history += S230
- `archive/sessions/session230_commit_p1_composite_q_multi_axis_fft.md` — this file
- `status/CLOSED_PATHS.md` — slot-4 row appended
- `status/SESSION_INSIGHTS.md` — S230 entry appended
- `RESEARCH_AGENDA.md` — Arc 8 slot 4 marked done, slot 5 next
- `CROSS_DOMAIN_TECHNIQUES.md` — CRT-based multi-axis Dirichlet
  character DFT row added
- `.run_state` → 398

## Session-end self-evaluation

1. **What did I produce that was not in the project before this session?**
   - First end-to-end Dirichlet-L zero finder for COMPOSITE q in the
     project, via CRT-based multi-axis FFT primitive over direct-product
     groups (Z/p₁Z)* × … × (Z/p_KZ)*.
   - Multi-axis FFT primitive verified to 1e-15 relative error at
     q ∈ {15, 35, 105, 1001}.
   - mpmath cross-validation at q ∈ {15, 35} to ~0.05 abs accuracy.
   - End-to-end π(x; q, a) timing data at composite q.
   - Empirical FFT vs DIRECT benchmark at composite q showing 1.75× at
     q=1001 (vs slot-3's 2.04× at q=1009).
   - Methodological win: caught and fixed χ(-1) parity formula bug
     for composite q via mpmath cross-check.

2. **What edges did my work compose or cite?**
   E1.5 (per-query bit-content), E2.1 (spectral/sieve, multi-axis DFT
   generalises slot-2 cyclic DFT to direct-product groups via CRT),
   E3.1 (CCM setup-vs-eval), E6.6 (Aggarwal orthogonal), E6.7 (DR
   sieve barrier orthogonal), S224 (Correlation Dichotomy template),
   S227, S228, S229.

3. **If my session produced only duplicate closures, why?**
   It didn't. Multi-axis FFT primitive is a novel construction in
   the project; CRT-based dual-group DFT was not previously used.
   The end-to-end zero finder for composite q is operationally new.

4. **What is the next-action for the next agent?**
   Slot 5 (theoretical wrap):
   - PRIORITY (a): cross-conductor amortisation argument (slot-3
     leftover priority-a). Likely structurally impossible by CRT
     independence across distinct conductors; rigorous lower bound
     would seal Thread 6 negatively across all four amortisation
     axes.
   - PRIORITY (c): Riemann-Siegel correction terms for Dirichlet L
     (Berry 1995 / Coffey 2003), deferred from slot 4.
   - Synthesis: position Thread 6 results within Threads 1-5 closures
     (S190, S202, S215, S224). Final result is most likely the
     fourth Thread 6 amortisation-axis closure: a-direction trivial,
     χ-direction constant lift, q-axis composite extension constant
     lift, Q-batched structurally impossible.

   See `slot4_composite_q_multi_axis_fft_results.md` falsifiers F4.1-F4.4
   and open list.

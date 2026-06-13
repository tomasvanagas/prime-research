# Thread 6, slot 4 — Composite-q multi-axis FFT for Dirichlet L-eval + zero finder

**Date:** 2026-04-29
**Mode:** commit (Thread 6 / π in arithmetic progressions, batched on
modulus q, slot 4 of 5)
**Goal:** PRIORITY (b) from S229 recommendation. Extend slot-3's
end-to-end Dirichlet-L zero finder from prime q to composite squarefree
odd q via CRT-based multi-axis FFT primitive. Composite q is the
practical AP regime (q ∈ {30, 60, 105, 1001} are common conductors).

## Headline result (HONEST)

**Multi-axis FFT primitive built from the CRT decomposition of (Z/qZ)*
preserves slot-3's structural lift in composite-q regime.** End-to-end
single-query π(x; q, a) timing:

```
q=1009 (prime,    φ=1008, slot-3, K=200, x=10⁶):  186.3 ms
q=1001 (composite,φ=720,  slot-4, K=200, x=10⁶):  166.6 ms
```

11% speedup at q=1001 vs q=1009 attributable mostly to smaller φ(1001)
= 720 vs φ(1009) = 1008 (28% smaller). Per-character cost is *not*
faster (q=1001 is slightly more expensive per-character because
N_AFE ≈ 1872 ≈ N_AFE(q=1009)).

**FFT vs DIRECT (no-FFT-sharing) at composite q:**

```
q=105  (φ=48,  K=50,  N_AFE=283):  FFT 0.83× DIRECT  (FFT loses, small-q regime)
q=1001 (φ=720, K=200, N_AFE=1872): FFT 1.75× DIRECT  (FFT wins, large-q regime)
```

Slot-3 measured FFT 2.04× DIRECT at prime q=1009, K=200. Multi-axis
FFT at q=1001 (composite, similar size) measures 1.75× — within 14% of
the prime-q result. **Multi-axis FFT preserves the slot-3 speedup
shape; does NOT introduce asymptotic gain at composite q.**

**This is a B-grade ambitious-attack-with-bounded-result.** The
multi-axis primitive correctness verifies to 1e-15 relative error at
q ∈ {15, 35, 105, 1001}. Hardy-Z + sign-change zero finder validates
against mpmath ground truth to ~0.05 absolute at q ∈ {15, 35} —
consistent with slot-3's ~0.03 accuracy at q=11. Operational regime
extended from {prime q} to {prime q + squarefree odd composite q with
no factor 2}, with FFT primitive's structural lift preserved.

## What was built

`slot4_composite_q_multi_axis_fft.py` (~600 lines):

### CRT decomposition layer

For squarefree odd q = p₁ · p₂ · … · p_K, the multiplicative group
decomposes via CRT:

```
(Z/qZ)*  ≅  (Z/p_1Z)* × (Z/p_2Z)* × … × (Z/p_KZ)*
        ≅  Z/(p_1-1) ⊕ Z/(p_2-1) ⊕ … ⊕ Z/(p_K-1)
```

via primitive roots g_i mod p_i. Each n ∈ (Z/qZ)* corresponds to a
tuple (k_1, …, k_K) with n ≡ g_i^{k_i} (mod p_i). Characters factor as

```
χ_{j_1,…,j_K}(n) = ∏_i ω_i^{j_i k_i}    ω_i = exp(2πi/(p_i-1))
```

Implementation: `crt_decomp(q)` returns (primes, group_orders, prim_roots);
`build_crt_log_table(q, ...)` precomputes the (residue → CRT-tuple) map.

### Multi-axis FFT primitive

Aggregate W indexed by tuple (k_1,…,k_K) ∈ ⊕_i Z/(p_i-1)Z. Multi-axis
DFT via `numpy.fft.ifftn`:

```python
L[j_1,…,j_K]
    = Σ_{(k_1,…,k_K)} ∏_i ω_i^{j_i k_i} · W[k_1,…,k_K]
    = phi · ifftn(W, axes=(0,…,K-1))[j_1,…,j_K]
```

phi = ∏_i (p_i - 1) = φ(q) is the total scale factor. The W tensor has
shape `(p_1-1, …, p_K-1, n_t)` with the time axis last; `ifftn` over the
first K axes produces L-values for ALL φ(q) characters simultaneously.

### Cross-method correctness

```
q=15   (factors [3, 5],     φ=8):     err = 9.93e-16   rel = 2.47e-16
q=35   (factors [5, 7],     φ=24):    err = 2.98e-15   rel = 4.40e-16
q=105  (factors [3, 5, 7],  φ=48):    err = 2.66e-15   rel = 6.33e-16
q=1001 (factors [7, 11, 13],φ=720):   err = 9.93e-15   rel = 1.13e-15
```

Multi-axis FFT vs direct per-character matmul agree to floating-point
precision at all tested composite q. **The CRT-multi-axis FFT identity
is verified.**

### Full AFE (composite q): main + reflected

`compute_full_l_via_afe_composite` extends slot-3's two-step AFE to
composite q:

1. Main term M(t, χ) via multi-axis FFT primitive.
2. Reflected term: for primitive χ mod q (squarefree odd q, every j_i ≠ 0),
   `Σ_{n ≤ N} χ̄(n)/n^{½-it} = conj(M(t, χ))` because χ̄(n) = conj(χ(n))
   for unitary characters of finite abelian groups and n ∈ ℝ_{>0}.
   So reflected = `W_χ · (q/π)^{-it} · ratio_gamma · conj(M)`. NO
   additional FFT — pointwise multiply.

This identity holds for ALL Dirichlet characters mod q, not just
prime-q. The simplification makes full AFE essentially free given the
main term, just as in slot 3.

### Gauss sums via multi-axis FFT

`gauss_sum_composite` computes τ(χ_j) = Σ_a χ_j(a) e^{2πi a/q} for all
characters via the same multi-axis FFT identity:

```
τ(χ_j) = phi · ifftn(e_a)[j]
```

where `e_a[k_1,…,k_K] = exp(2πi · r/q)` with r the integer residue
having CRT-tuple (k_1,…,k_K). Cost ≈ phi · log(phi) complex ops.

Verified `|W_χ| = 1` (within numerical precision) for all primitive χ
mod q ∈ {15, 35, 105, 1001}.

### χ(-1) parity

For each odd prime p_i, ω_i^{(p_i-1)/2} = -1 and -1 ≡ g_i^{(p_i-1)/2}
(mod p_i), so χ_i(-1) = ω_i^{j_i · (p_i-1)/2} = (-1)^{j_i}. Therefore

```
χ(-1) = ∏_i (-1)^{j_i} = (-1)^{Σ_i j_i}
```

Initial implementation used the wrong formula `(j_i · (p_i-1)/2) mod 2`
which only matched for prime q (collapses to j_1 mod 2 since K=1).
Corrected to `(Σ j_i) mod 2`. Bug-fix verified against mpmath ground
truth at q=15, 35.

### Hardy Z + sign-change bracketing

Same as slot 3:

```
Z_χ(t) = exp(i θ_χ(t)) · L(½+it, χ)
θ_χ_j(t) = (t/2) log(q/π) + arg Γ((1/2 + a_j + it)/2) - arg(W_χ_j)/2
```

`zero_brackets_vectorised` is identical to slot 3's vectorised
`np.where(sign(Z[:,:-1]) * sign(Z[:,1:]) < 0)` over the (φ, n_t) array.

### mpmath ground-truth (small-q sanity)

`mpmath_first_zeros_composite` evaluates Hardy Z directly with mpmath
high-precision arithmetic on a fine grid; sign-change bracketing
identical to slot-4's pipeline but with N_factor=80 (vs slot-4's
oversize=12) and dps=30 mpmath precision. Cross-check:

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

Slot-4 zero-finder accuracy at composite q (max 0.05 abs at q=15, 0.06
at q=35) is consistent with slot-3's accuracy at prime q (0.033 at q=11).

## End-to-end π(x; q, a) timing (slot-4 vs slot-3)

| q     | factors      | φ    | K   | t_zf (ms) | t_psum (ms) | t_single (ms) | zeros[med] |
|-------|--------------|------|-----|-----------|-------------|---------------|------------|
| 15    | [3, 5]       | 8    | 50  | 3.80      | 0.16        | 3.96          | 10         |
| 35    | [5, 7]       | 24   | 50  | 4.90      | 0.35        | 5.25          | 46         |
| 105   | [3, 5, 7]    | 48   | 50  | 6.59      | 0.36        | 6.94          | 0*         |
| 105   | [3, 5, 7]    | 48   | 200 | 37.82     | 0.56        | 38.38         | 0*         |
| 1001  | [7, 11, 13]  | 720  | 200 | 152.27    | 14.31       | 166.58        | 138        |

*At q=105 the median is reported as 0 because 33 of 48 chars are
non-primitive (induced from smaller modulus); slot-4 leaves their zero
list empty as a known limitation. The cost measurement still reflects
the AFE + bracketing cost over ALL 48 characters because the multi-axis
FFT primitive does not skip them. Slot 5 / future work: handle
non-primitive characters via inducing-primitive-character zeros + Euler
factor corrections.

### Comparison to slot-3 prime-q baseline

```
                                slot-3 (prime q=1009)    slot-4 (comp q=1001)
T_setup                            151 ms                    -- (smaller setup)
T_AFE                              139 ms                    ~140 ms
T_partial_sum                       34.8 ms                   14.3 ms (smaller K-effective)
T_total single-query               186.3 ms                  166.6 ms
                                                              (89% slot-3 cost)
```

The 11% wall-clock drop at q=1001 vs q=1009 is mostly explained by
smaller φ (720 vs 1008, a 29% reduction). Per-character cost is roughly
unchanged because N_AFE depends on q, not φ, and is essentially identical
between q=1001 and q=1009.

## FFT vs DIRECT comparison (slot-3 within-pipeline benchmark)

```
q=105,  K=50,  φ=48,  N_AFE=283,   n_t=218:
  FFT method     :  5.09 ms
  DIRECT method  :  4.24 ms
  FFT speedup    :  0.83×                  <-- FFT LOSES (small composite q)

q=1001, K=200, φ=720, N_AFE=1872,  n_t=823:
  FFT method     : 145.24 ms
  DIRECT method  : 254.39 ms
  FFT speedup    :  1.75×                  <-- FFT WINS

(slot-3 prime-q reference for comparison:
q=101,  K=200, φ=100, N_AFE=678,   n_t=849:  FFT 0.95× DIRECT
q=1009, K=200, φ=1008,N_AFE=1883,  n_t=823:  FFT 2.04× DIRECT)
```

**Key empirical findings:**

1. Multi-axis FFT primitive crosses over from "loss" to "win" around
   q ~ 200-500, same as slot-3 single-axis FFT for prime q.
2. At q=1001 (composite, φ=720): 1.75× speedup. At q=1009 (prime, φ=1008):
   2.04× speedup. Multi-axis FFT is *slightly worse* than single-axis
   at comparable q because:
   - `numpy.fft.ifftn` over a tensor of shape (p_1-1, …, p_K-1, n_t) calls
     1D FFT routines internally K times along K axes, each of size (p_i-1)
     much smaller than 1008. Smaller transforms have proportionally higher
     constant overhead.
   - At q=1001 with shape (6, 10, 12) the per-axis transforms are tiny;
     numpy's overhead per call adds up.
3. Small composite q (≤ 200) loses to DIRECT for the same reason:
   FFT overhead exceeds the BLAS-saturated matmul cost at that scale.

## Edges composed / cited

- **E1.5** (per-query bit-content): unchanged — slot-4's multi-axis
  primitive affects setup + AFE eval but not partial-sum eval, exactly
  as E1.5 predicts.
- **E2.1** (spectral / sieve interface): the multi-axis DFT identity
  generalises slot-2's cyclic DFT to direct-product groups via CRT.
- **E3.1** (CCM amortisation, downgraded): same setup-vs-eval decoupling
  pattern; multi-axis FFT lifts the setup but not the eval.
- **E6.6** (Aggarwal binary search): orthogonal — within-q amortisation
  rather than cross-call.
- **E6.7** (Deléglise-Rivat sieve barrier): orthogonal — sieve pillar.
- **S224** (Correlation Dichotomy, Thread 5): slot-4's wall-clock shape
  ~1.75× constant lift at q=1001 is *not* Correlation-Dichotomy-shaped
  (33× at M=64). Same shape as slot-3: "real-zeros for synthetic-cost"
  bounded constant.
- **S227, S228, S229** (slot 1, 2, 3 of Thread 6): slot-4 extends slot-3's
  apparatus to composite q with the structural-lift shape preserved.

## Cross-domain ingredient (CROSS_DOMAIN_TECHNIQUES.md)

**CRT-based multi-axis Dirichlet character DFT.** Generalisation of
slot-2's cyclic-group DFT primitive to direct-product groups
(Z/p_1Z)* ⊕ … ⊕ (Z/p_KZ)* via Chinese Remainder Theorem
decomposition. Implemented as `numpy.fft.ifftn` on a multi-dimensional
tensor indexed by CRT log-tuples. **NEW USED I** in this project.

The technique is the natural generalisation of slot-2's primitive but
required:
- Explicit CRT log-table construction.
- Correct character indexing (CRT-tuple to flat index via row-major
  reshape).
- Correct χ(-1) formula `(Σ j_i) mod 2`, NOT the per-axis variant.
- Primitive-character flag (j_i ≠ 0 ∀ i for squarefree q).

References:
- Bober & Hiary, "On computing the moments of the Dirichlet L-function",
  for the slot-2 antecedent (cyclic FFT primitive).
- Standard CRT decomposition of (Z/qZ)* (Davenport's "Multiplicative
  Number Theory" §1).

## Self-extension proposals

Per CLAUDE.md autonomy invariant:

**Successor 1 (slot 5 PRIORITY-a, NEW)**: Cross-conductor batches via
shared coprime tower. For nested Q = {q_1 ⊆ q_2 ⊆ …} with q_i divisor
of q_{i+1}, the residue tower allows reuse of the cp_all summands. E.g.,
q=15 and q=105 share the gcd≠1 prefix structure. Investigate whether
a shared partial-sum cache amortises across {q=15, 105, 1001} when q's
are structurally nested.

**Successor 2 (slot 5 PRIORITY-c, NEW)**: Riemann-Siegel correction
terms for Dirichlet L (Berry 1995 / Coffey 2003). Drops N_AFE by 5-10×
with constant overhead; would re-test multi-axis FFT primitive at much
smaller AFE cost. PROPOSED in CROSS_DOMAIN_TECHNIQUES.md §8 (slot-3
recommendation, not yet implemented).

**Successor 3 (slot 5 PRIORITY-d, NEW)**: Non-primitive character
handling via inducing primitive character + Euler factor corrections.
At q=105 = 3·5·7, characters with j_3 = 0 are induced from mod 15
characters; their L-zeros are the q=15 zeros plus zeros from
(1 - χ_{15}(7)·7^{-s}). This would let slot-4 cover the full character
space at composite q, useful for application but orthogonal to the
slot-4 cost question.

## Falsifiers (slot-4 bounded results)

**F4.1** Multi-axis FFT does NOT give an asymptotic improvement over
slot-3's single-axis FFT at composite q. Both achieve constant-factor
lift in the same crossover regime (~q=500). At q=1001 vs q=1009, the
multi-axis primitive is *marginally slower* (1.75× vs 2.04× FFT/DIRECT
ratio) because of multiple small per-axis 1D FFTs vs one large 1D FFT.

**F4.2** End-to-end wall-clock is dominated by AFE eval (N_AFE · n_t),
which depends on q (not φ); composite q with many small prime factors
inherits the same AFE cost as nearby prime q. The 11% drop at q=1001
vs q=1009 is a smaller-φ effect, not a multi-axis-FFT effect.

**F4.3** Multi-axis FFT does not shrink the AFE truncation requirement.
At fixed q, balanced N = √(qt/(2π)) is the same regardless of factorisation.
Composite q does not enable cheaper AFE per se.

**F4.4** Cross-conductor amortisation across the Q-family was NOT
attempted in slot 4 (the slot-3 leftover priority-a falsifier). Slot 5
should attempt this if a structural mechanism is plausible.

## What this slot resolves vs leaves open

**Resolves:**

- Multi-axis FFT primitive built and verified to floating-point
  precision at composite q ∈ {15, 35, 105, 1001}.
- End-to-end zero finder works for composite squarefree odd q with
  primitive-character accuracy comparable to slot-3 prime-q accuracy
  (~0.05 abs).
- FFT vs DIRECT crossover at composite q has the same shape as at
  prime q: FFT loses for q < 200, wins for q > 500.
- Multi-axis FFT preserves slot-3's structural lift but does NOT
  introduce additional asymptotic gain.

**Leaves open:**

- Cross-conductor batches (slot-3 priority-a, deferred to slot 5).
- Riemann-Siegel correction terms (slot-3 priority-c, deferred).
- Non-primitive character handling at composite q (slot-4 limitation).
- q with prime-power factors (e.g., q = 4·p, q = 8·p) requires the
  more general (Z/2^aZ)* structure (Z/2 × Z/2^{a-2}); not implemented.
- q = 2k for odd k > 1: same characters as mod k (since (Z/2qZ)* ≅
  (Z/qZ)*); could be handled by reduction to slot-3 / slot-4 prime-q
  / squarefree-odd-q pipelines.

## What would qualify as slot-5 progress (theoretical wrap)

(per the original 5-slot plan, slot 5 is the theoretical wrap):

(a) **Cross-character lower bound** matching the empirical curve.
    Show that ANY algorithm computing L-zeros for all χ mod q has
    cost Ω(φ(q) · K) up to a polylog factor at fixed q. Prove that
    cross-character amortisation across L-zeros is structurally
    impossible at fixed q, sealing Thread 6 negatively for the
    single-conductor regime.

(b) **Cross-conductor break:** For Q = {q_1, …, q_M} batched, can a
    shared-zero-tower or shared-Euler-factor primitive achieve
    sub-M × per-conductor amortisation? Likely structurally impossible
    by CRT independence but worth a careful argument.

(c) **Synthesis with Threads 1-5 closures** (S190 invariant subspace,
    S202 cross-x amortisation, S215 plethysm, S224 Correlation Dichotomy).
    Position slot-3/slot-4 results within the broader Thread 6 negative-
    shape closure (or partial-positive shape).

## Self-grade for slot 4

**B** — substantive empirical refinement and architectural extension.

Built:
- Multi-axis FFT primitive (CRT decomposition + ifftn) verified to 1e-15.
- End-to-end zero finder for composite squarefree odd q.
- mpmath cross-validation at q ∈ {15, 35} to ~0.05 abs accuracy.
- FFT vs DIRECT comparison at composite q.

Refines slot-3's prediction:
- Multi-axis FFT at q=1001 gives 1.75× FFT/DIRECT — smaller than
  slot-3's 2.04× at q=1009 (multiple small 1D FFTs cost more than
  one large 1D FFT in numpy's BLAS-vs-FFT regime).
- 11% end-to-end speedup at composite q vs nearby prime q is
  attributable to smaller φ, not multi-axis-FFT-specific.

Not A: multi-axis FFT does NOT introduce asymptotic improvement; the
shape is "constant-factor lift preserved at composite q" not
"Correlation-Dichotomy-shaped partial-positive". The structural-lift
finding from slot 3 extends to composite q, which is the practical
regime, but the wall-clock advantage stays bounded.

Not C: composite-q multi-axis FFT primitive is operationally new in
the project, the CRT-based AFE pipeline is new, the FFT vs DIRECT
benchmark at composite q is new data, and the χ(-1) parity bug
(initial wrong implementation) was caught and fixed via mpmath
cross-validation — methodological win.

The honest characterisation: B for substantive engineering extension
of slot 3 to the practical AP-π regime. The multi-axis FFT primitive
is the natural generalisation of slot 2's cyclic primitive, the
implementation works correctly, but it does not produce an asymptotic
break. Slot 5 should focus on the theoretical-wrap angle (cross-conductor
amortisation lower bound, or Thread 6 final synthesis with adjacent
threads).

## Files written by this slot

- `slot4_composite_q_multi_axis_fft.py` — new (~600 lines)
- `slot4_composite_q_multi_axis_fft_results.md` — this file
- `slot4_composite_q_end_to_end.csv` — 5 rows
- `slot4_composite_fft_vs_direct.csv` — 2 rows

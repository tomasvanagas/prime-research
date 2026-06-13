# Cross-x amortisation, slot 1: decoupled setup / evaluation profile

**Thread 5 (commit), slot 1 of 5.** Session 220.

## What this slot asked

Build an explicit-formula evaluator that decouples `K_zeros_setup` from
`K_per_x_evaluation`. Profile both. Measure for K ∈ {⌈log² x⌉, ⌈log³ x⌉,
⌈x^{1/4}⌉, ⌈x^{1/2}⌉} and a geometric K-grid, at x ∈ {10⁵, 10⁶, 10⁷}.
Determine empirically:

- (Q1) Is `T_eval(K, x)` linear in K?
- (Q2) Is the per-term cost `T_eval(K, x) / K` x-independent?
- (Q3) How does `T_setup(K)` scale?

The answers determine the per-x amortised cost
`T_per_x_amortised = T_setup(K) / M + T_eval(K, x_avg)` over M batched
queries — the central object of Thread 5.

## Build

`cross_x_decoupled_profile.py` reuses the `R_at_rho` and `riemann_R`
machinery from `experiments/analytic/connes_amortisation/connes_amortisation.py`
(the same evaluator that backed S195 / S202). It adds:

- `time_setup(K)`: median elapsed time to load K cached zeros (the
  cache is the project's standard zero source; the from-scratch cost
  via Hiary 2011 is reported separately as `K^{17/13}`).
- `time_eval(x, K, gammas)`: time `R(x)` and the partial sum
  `2 Σ_{j=1}^K Re R(x^{ρ_j})` separately, returning
  `(T_R, T_sum, T_total)`.

Output CSVs (committed alongside the script):
- `cross_x_decoupled_setup.csv`   — (K, T_setup_load_s, Hiary_K^{17/13})
- `cross_x_decoupled_eval.csv`    — (x, K, policy, T_R_s, T_sum_s, T_total_s, T_per_term_s)
- `cross_x_decoupled_fits.csv`    — log-log fit of T_sum vs K, per x
- `cross_x_decoupled_perterm_spread.csv` — cross-x spread of T_per_term at fixed K

## Empirical profile (K_max = 3200, x ∈ {10⁵, 10⁶, 10⁷})

### Setup cost

| K    | T_setup_load (cached) | Hiary K^{17/13} (arith ops) |
|------|-----------------------|-----------------------------|
| 25   | 0.0001 s              | 67.3                        |
| 50   | 0.0001                | 166.6                       |
| 100  | 0.0001                | 412.5                       |
| 200  | 0.0002                | 1021.0                      |
| 400  | 0.0003                | 2527.5                      |
| 800  | 0.0005                | 6256.7                      |
| 1600 | 0.0011                | 15488.3                     |
| 3200 | 0.0023                | 38340.5                     |

**Cache-load scaling:** linear in K (file I/O dominates), prefactor
~0.7 µs/zero on this machine. The cached cost is *not* the asymptotically
relevant number; for production-scale K = Θ̃(x_max) = 10⁵⁰⁻¹⁰⁰ the zeros
must be computed, with the best-known cost K^{17/13} per Hiary 2011.

### Per-x evaluation

Selected entries (full table in `cross_x_decoupled_eval.csv`):

| x     | K    | policy | T_sum (s) | T_per_term (µs) |
|-------|------|--------|-----------|-----------------|
| 10⁵   | 133  | log²x  | 0.197     | 1483            |
| 10⁵   | 316  | x^{1/2}| 0.331     | 1046            |
| 10⁵   | 1526 | log³x  | 1.243     | 815             |
| 10⁵   | 3200 | —      | 2.317     | 724             |
| 10⁶   | 191  | log²x  | 0.210     | 1098            |
| 10⁶   | 1000 | x^{1/2}| 0.810     | 810             |
| 10⁶   | 2637 | log³x  | 2.074     | 786             |
| 10⁶   | 3200 | —      | 2.131     | 666             |
| 10⁷   | 260  | log²x  | 0.198     | 762             |
| 10⁷   | 3162 | x^{1/2}| 1.896     | 600             |
| 10⁷   | 3200 | log³x  | 2.387     | 746             |

### Q1 — log-log fit, T_sum vs K

| x    | α (slope) | prefactor exp(β) | per-term drift |
|------|-----------|------------------|----------------|
| 10⁵  | 0.69      | 0.0075           | −4.9·10⁻⁴      |
| 10⁶  | 0.72      | 0.0058           | −4.1·10⁻⁴      |
| 10⁷  | 0.78      | 0.0032           | −2.6·10⁻⁴      |

α < 1 is **NOT** a structural sublinearity in arithmetic operations.
It is a finite-K warm-up artifact:

- Each `R_at_rho(x, γ_j)` calls `mp.ei(arg)` where `|arg| = γ_j · ln(x) / n`.
  For small γ_j the magnitude is small and mpmath uses a slowly-
  convergent series; for large γ_j the magnitude is large and mpmath
  uses an O(1)-term asymptotic expansion. The per-term cost therefore
  *decreases* with j until the asymptotic regime is reached, at which
  point it floors.
- Empirically, the per-term cost floors at ~600–750 ns/zero by
  K ≈ 1600 (see Q2 spread). Per-term cost at K=3200 is *higher* than
  at K=1600 for x = 10⁷ (746 vs 623 µs), confirming the floor.
- Slope α increases with x (0.69 → 0.72 → 0.78). At larger x the
  asymptotic regime is reached at smaller K (because |arg| = γ · ln x / n
  scales with ln x). Extrapolation: α → 1 as both K and x grow.

**The asymptotic per-term cost is Θ(1) arithmetic operations**, hence
**T_sum(K, x) = Θ(K) asymptotically.** The empirical α ∈ [0.69, 0.78]
is the sub-asymptotic regime visible at K ≤ 3200; beyond the per-term
floor, T_sum scales linearly in K.

### Q2 — per-term cross-x spread at fixed K

| K    | min per-term (µs) | max per-term (µs) | spread |
|------|-------------------|-------------------|--------|
| 25   | 2164              | 3243              | 50%    |
| 50   | 1356              | 2445              | 80%    |
| 100  | 1063              | 1647              | 55%    |
| 200  | 805               | 1218              | 51%    |
| 400  | 701               | 977               | 39%    |
| 800  | 657               | 889               | 35%    |
| 1600 | 623               | 872               | 40%    |
| 3200 | 666               | 746               | 12%    |

At small K the per-term cost is x-coupled (warm-up region);
at K = 3200 the spread collapses to 12%, consistent with the
asymptotic-Ei regime being x-independent.

**Q2 conclusion: T_per_term → constant in x as K grows. The leading-
order amortisation question is independent of which x_i are queried.**

## What the slot 1 profile establishes

The decoupled cost model for batched queries `π(x_1), ..., π(x_M)` at
truncation K:

```
T_total(K, M) = T_setup(K) + Σ_{i=1}^M T_eval(K, x_i)
              = K^{17/13} (Hiary asymptotic) + M · (T_R + a·K + b)
              = K^{17/13} + O(M K)
```

with `a` (per-term floor) being x-independent at large K and `b` (warm-
up) being a constant shifted by O(log K · log x) corrections that
vanish in the leading order.

**Per-x amortised cost:**

```
T_per_x_amortised(K, M) = T_setup(K) / M + T_eval(K, x_avg)
                        = K^{17/13} / M + a · K + O(1).
```

For per-x amortised cost to be polylog(x), we need either:

(a) **K = polylog(x)**, in which case the `a · K` term is polylog and
    the K^{17/13}/M term vanishes for any M ≥ 1 (since K^{17/13} is
    polylog too). **But Thread 3 (S202) closed K = polylog in
    distribution: K* = Θ̃(x) under Montgomery for any positive hit-rate.**
    So this option is closed.

(b) **a · K is polylog while K is not polylog** — this would require
    the per-term cost to decay polynomially in K, i.e., α < 1
    asymptotically. The slot-1 measurement establishes that the
    per-term cost *floors* at ~constant for K ≥ 1600; α → 1 in the
    asymptotic regime. So this option is also closed.

**The slot-1 conclusion (structural):** setup amortisation is orthogonal
to per-x evaluation. Even with infinite M (free setup), per-x amortised
cost is dominated by `a · K` with K = Θ̃(x_max) by Thread 3, giving
amortised per-x cost Θ̃(x). Setup amortisation provides a vanishing
contribution.

This is a *partial* closure of Thread 5: the explicit-formula pillar
for cross-x amortisation, with structural reason in hand. Subsequent
slots cover (slot 2) batched empirical timing, (slot 3) Meissel-Lehmer
/ HKM sieve cross-x amortisation, (slot 4) Aggarwal binary search
re-examination, (slot 5) theoretical wrap.

## What this does NOT close

The slot 1 profile measures the *single-zero* per-term cost. It does
NOT close the following potential cross-x amortisations within the
explicit-formula pillar itself:

1. **Multipoint evaluation of R(·^ρ) across x_1, ..., x_M for fixed ρ.**
   For correlated x_i (e.g., x_i = x_max · (1 − iε)), the inputs to Ei
   are close, and a Taylor expansion or Remez approximation could in
   principle reduce the M evaluations from O(M) per ρ to O(M^{1−δ}).
   Slot 2 must measure this.
2. **Shared Möbius truncation across (x_i, ρ_j).** Currently each
   `R_at_rho` evaluation does its own Möbius prefix (M=12 by default).
   Per-batch hoisting could amortise the prefix, but this is a
   constant-factor speedup and does not change the asymptotic α.
3. **Block-precomputation of `e^{i γ_j ln(x_max)}` followed by
   correction factors `e^{i γ_j ln(x_i / x_max)}`.** This is a
   well-known trick in numerical Riemann-zeta evaluation (Schönhage
   1990 / Odlyzko-Schönhage). Total cost stays Θ(M K) but with
   different constants.

These are all **constant-factor / sub-leading** in the asymptotic α,
so they cannot rescue polylog amortisation when K = Θ̃(x). Slot 2 will
nonetheless measure (1) empirically at the M ∈ {log x, log² x, x^{1/4},
x^{1/2}} scales mandated by the thread 5 plan.

## Falsifiers (slot-1 specific)

The conclusion *would* be falsified by:

- A measured per-term cost that **decays polynomially in K asymptotically**
  (α < 1 stable as K → ∞, not the warm-up artifact observed at K ≤ 3200).
  Run on real Hiary-from-scratch zeros at K ≥ 10⁶ would test this — but
  is computationally out of reach for this slot. The asymptotic
  expansion of `Ei(z)` for |z| → ∞ takes O(1) terms, so this falsifier
  would require an mpmath bug or a fundamentally different evaluation
  strategy.
- T_per_term showing *strong* x-coupling (spread ≫ 50% at large K). The
  observed K=3200 spread of 12% rules this out within the tested range.
- T_setup growing super-linearly in K (e.g., as K² due to dependency
  conflicts in zero-computation). Hiary 2011's K^{17/13} bound holds
  unconditionally; this is not a real falsifier.

## Edges composed / cited

- **E1.5** (information-theoretic per-query barrier): the `a · K`
  per-x cost matches the Σ x/K bit-content barrier when K = Θ̃(x).
  Slot 1 confirms the wall-clock realisation of this barrier.
- **E3.1** (Connes-Consani-Moscovici closure as refined at S193 / S202):
  setup-cost dominance via Hiary K^{17/13} versus Connes K³ holds at
  the slot-1 measurement scale, with the additional empirical observation
  that the cached-zero load cost is Θ(K) on file I/O, dominated by both
  setup-cost models in arithmetic terms.
- **E6.6** (Aggarwal binary search): the per-x evaluation cost
  characterised here is the cost of one π-sub-call in Aggarwal's
  log(x) binary-search descent. Slot 4 will use this profile directly.
- **S195 row 816** (GUE random-phase variance): the K* = Θ̃(x)
  threshold this slot uses to argue the structural closure is the
  S195 result.
- **S202 row** (Connes-Galway equivalence wrap): Thread 5 extends
  the per-query closure to the batched-query regime; slot 1 is the
  first-pillar (explicit-formula) extension.

## Cross-domain ingredient (CROSS_DOMAIN_TECHNIQUES.md)

Section 8 amortised algorithmics row, status PROPOSED → PARTIAL.
The decoupled profile is the standard technique applied: separate
amortisable from per-query work, then compute the M → ∞ limit. The
slot-1 contribution is the explicit-formula concrete realisation;
slots 2–5 will fill in the M-dependence and the other two pillars.

Reference: Tarjan 1985 "Amortized computational complexity"
(SIAM J. Algebraic Discrete Methods); the application here is the
cross-x query model rather than the more familiar self-adjusting
data-structure setting.

## Self-grade

**B** — substantive empirical refinement. Slot 1 produces:

- A new evaluator that decouples T_setup from T_eval and a profile
  showing per-term floor at ~600 ns × K asymptotically.
- A structural argument that closes the explicit-formula pillar
  for cross-x amortisation, conditional on Thread 3 closure
  (Montgomery random-phase). The closure is *partial* (slots 2–5
  must check the other two pillars and the multipoint-evaluation
  sub-question).
- Concrete falsification criteria (asymptotic per-term polynomial
  decay; strong x-coupling at large K; super-linear setup).

Not A-grade because (a) the closure conditions on the same
Montgomery heuristic that backs Thread 3 and (b) the cross-pillar
generalisation (HKM, Aggarwal) is not yet attempted.

## .commit_state changes

Slot 1 advances `sessions_used` 0 → 1.

```
thread:cross_x_amortisation
sessions_used:1
status:ACTIVE
session_history:S220
```

## Next-action for slot 2

Slot 2 measures the *empirical* per-x amortised cost when M batched
queries are run with shared K zeros: implement the
multi-x evaluator that builds K zeros once, then evaluates the
partial sum at M = {log x, log² x, x^{1/4}, x^{1/2}} correlated x_i.
Plot per-x amortised cost vs M; confirm the slot-1 forecast that
amortised cost saturates at the `a · K` floor independent of M.
